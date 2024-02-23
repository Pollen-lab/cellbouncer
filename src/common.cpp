#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <mixtureDist/functions.h>
#include <htswrapper/bc.h>

/**
 * Contains functions used by more than one program in this
 * repository.
 */

using std::cout;
using std::endl;
using namespace std;

/**
 * Parse a file output by demux_mt or demux_vcf, and store barcodes mapped
 * to individual IDs.
 */
void parse_barcode_map(string& fn, 
    map<unsigned long, string>& bc2hap,
    set<string>& barcode_groups,
    double llr_cutoff,
    bool keep_doublets){

    ifstream infile(fn.c_str());
    
    string bc_str;
    string hap_str;
    char singdoub;
    double llr;

    while (infile >> bc_str >> hap_str >> singdoub >> llr){
        if ((keep_doublets || singdoub == 'S') && llr >= llr_cutoff){
            // Trim any trailing stuff off of the barcode, so it can
            // be interpreted as a bitset
            unsigned long bc_hashed = bc_ul(bc_str); 
            bc2hap.insert(make_pair(bc_hashed, hap_str));
            barcode_groups.insert(hap_str);
        }
    }
}

/**
 * For doublet identification (in a data structure), convert
 * a combination of haplotype indices i and j into a single index.
 */
short hap_comb_to_idx(short i, short j, short nhaps){
    if (i >= nhaps || j >= nhaps){
        return -1;
    }
    if (i > j){
        // Swap.
        short tmp = j;
        j = i;
        i = tmp;
    }
    // The overall index in the data structure (unique idx for combination)
    short combined_idx = nhaps;
    // The "i" haplotype in the combination
    short first_idx = 0;
    // The "j" haplotype in the combination
    short second_idx = 1;
    short k = 1;
    while (first_idx < i){
        combined_idx += (nhaps - k);
        first_idx++;
        second_idx++;
        k++;
        if (k >= nhaps){
            return -1;
        }
    }
    // first_idx should now == i. Find j
    while (second_idx < j){
        second_idx++;
        combined_idx++;
    }
    return combined_idx;
}

/**
 * Undo the above combination - convert a combined haplotype idx (for doublets)
 * into a <i, j> pair denoting the two indices that went into creating it.
 */
pair<short, short> idx_to_hap_comb(short idx, short nhaps){
    if (idx < nhaps){
        return make_pair(-1, -1);
    }
    short combined_idx = nhaps;
    short i = 0;
    short j = 1;
    short k = 1;
    while (combined_idx < idx){
        if (idx > combined_idx && idx < combined_idx + (nhaps - k)){
            // Increment j till we find it.
            j += (idx - combined_idx);
            return make_pair(i, j);
        }
        combined_idx += (nhaps - k);
        j++;
        i++;
        k++;
        if (k >= nhaps){
            return make_pair(-1, -1);
        }       
    }
    if (combined_idx == idx){
        return make_pair(i, j);
    }
    return make_pair(-1, -1);
}

/**
 * Given a numeric index (single or doublet combination) and a vector
 * of (string) sample names, returns the name of the given sample
 * or sample combination.
 *
 * Ensures that names of doublet combinations are always given in 
 * alphabetic order.
 */
string idx2name(int x, vector<string>& samples){
    string indv_name;
    if (x < samples.size()){
        indv_name = samples[x];
    }
    else{
        pair<short, short> hc = idx_to_hap_comb(x, samples.size());
        if (samples[hc.first] < samples[hc.second]){
            indv_name = samples[hc.first] + "+" + samples[hc.second];
        }
        else{
            indv_name = samples[hc.second] + "+" + samples[hc.first];
        }
    }
    return indv_name;
}

/**
 * Initialize a pre-declared distance matrix by setitng all elements
 * that will be used and accessed to 0 and all others to -1.
 */
void init_distmat(vector<vector<float> >& dist_mat, int dim){
    for (int i = 0; i < dim; ++i){
        vector<float> row;
        for (int j = 0; j <= i; ++j){
            row.push_back(-1);
        }
        for (int j = i +1; j < dim; ++j){
            row.push_back(0);
        }
        dist_mat.push_back(row);
    }
}

/**
 * Print distance matrix to stderr.
 */
void print_distmat(vector<vector<float> >& dist_mat){
    for (int i = 0; i < dist_mat.size(); ++i){
        for (int j = i + 1; j < dist_mat.size(); ++j){
            if (j > i + 1){
                fprintf(stdout, " ");
            }
            fprintf(stdout, "%0.2f", dist_mat[i][j]);
        }
        fprintf(stdout, "\n");
    }
}

void print_distmat_square(vector<vector<float> >& dist_mat){
    for (int i = 0; i < dist_mat.size(); ++i){
        for (int j = 0; j < i; ++j){
            if (j != 0){
                fprintf(stdout, "\t");
            }
            fprintf(stdout, "%0.2f", dist_mat[j][i]);
        }
        fprintf(stdout, "\t0.00\t");
        for (int j = i + 1; j < dist_mat.size(); ++j){
            fprintf(stdout, "\t%0.2f", dist_mat[i][j]);
        }
        fprintf(stdout, "\n");
    }
}

/**
 * Given a table of log likelihood ratios between each possible pair of 
 * identities for a cell, iteratively eliminates the least likely possibility
 * (the individual deemed less likely by the highest-magnitude LLR in the table)
 * until only two individuals are left. The more likely of those individuals
 * is chosen as the true individual, and the LLR between the two final individuals
 * is chosen as the LLR of assignment.
 */
int collapse_llrs(map<int, map<int, double> >& llrs, double& llr_final){
    bool done = false;
    int idx_elim = -1;
    while (!done){
        int min_idx = -1;
        int max_idx = -1;
        double min_llr;
        int ncomps = 0;
        for (map<int, map<int, double> >::iterator llr = llrs.begin(); llr != llrs.end(); ){
            if (idx_elim != -1 && llr->first == idx_elim){
                llrs.erase(llr++);
            }
            else{
                for (map<int, double>::iterator llr2 = llr->second.begin(); 
                    llr2 != llr->second.end(); ){
                    if (idx_elim != -1 && llr2->first == idx_elim){
                        llr->second.erase(llr2++);   
                    } 
                    else{
                        if (min_idx == -1){
                            if (llr2->second < 0){
                                min_idx = llr->first;
                                max_idx = llr2->first;
                                min_llr = llr2->second;
                            }
                            else{
                                min_idx = llr2->first;
                                max_idx = llr->first;
                                min_llr = -llr2->second;
                            }
                        }
                        else if (abs(llr2->second) > abs(min_llr)){
                            if (llr2->second < 0){
                                min_idx = llr->first;
                                max_idx = llr2->first;
                                min_llr = llr2->second;
                            }
                            else{
                                min_idx = llr2->first;
                                max_idx = llr->first;
                                min_llr = -llr2->second;
                            }
                        }
                        ++ncomps;
                        ++llr2;
                    }
                }
                ++llr;
            }
        }

        if (ncomps == 0 || min_idx == -1){
            llr_final = 0;
            return -1;
        }
        else if (ncomps == 1){
            llr_final = -min_llr;
            return max_idx;
        }
        else{
            idx_elim = min_idx;   
        }
    }
    // Nothing found
    llr_final = 0.0;
    return -1;
}

/**
 * Given counts of identifications in a sample, compares the 
 * counts of different doublet combinations to expectations,
 * baed on the frequencies of each individual in the sample.
 *
 * Returns a chi-squared goodness of fit p-value.
 *
 * Low p-values indicate possibly incorrect doublet 
 * identifications, which suggests demultiplexing was
 * inaccurate.
 */
double doublet_chisq(map<int, int>& idcounts, int n_samples){
    
    if (n_samples <= 2){
        // Can't do test with up to only one doublet type
        return -1.0;
    }

    int tot_single = 0;
    int tot_double = 0;
    map<int, int> singles;
    map<int, int> doubles;
    for (map<int, int>::iterator ic = idcounts.begin(); ic != idcounts.end();
        ++ic){
        if (ic->first < n_samples){
            tot_single += ic->second;
            singles.insert(make_pair(ic->first, ic->second));
        }
        else{
            tot_double += ic->second;
            doubles.insert(make_pair(ic->first, ic->second));
        }
    }
     
    // If no doublets, can't do anything
    if (tot_double == 0){
        return 1.0;
    }

    // Get frequency of each single combination
    map<int, double> singfreq;
    for (map<int, int>::iterator s = singles.begin(); s != singles.end(); ++s){
        singfreq.insert(make_pair(s->first, (double)s->second/(double)tot_single));
    }
    
    // Check for missing doublet combinations, and store a count of 0 for each 
    for (int i = 0; i < n_samples-1; ++i){
        for (int j = i + 1; j < n_samples; ++j){
            int k = hap_comb_to_idx(i, j, n_samples);
            if (doubles.count(k) == 0){
                doubles.insert(make_pair(k, 0));
            }
        }
    }

    // Get expectation of each doublet combination
    map<int, double> doubfreq;
    double doubfreq_tot = 0.0; // Need to re-scale probs since not including self+self doublets
    for (map<int, int>::iterator d = doubles.begin(); d != doubles.end(); ++d){
        pair<int, int> combo = idx_to_hap_comb(d->first, n_samples);
        double expected = singfreq[combo.first] * singfreq[combo.second];
        doubfreq_tot += expected;
        doubfreq.insert(make_pair(d->first, expected));
    }
    for (map<int, double>::iterator df = doubfreq.begin(); df != doubfreq.end(); ++df){
        df->second /= doubfreq_tot;
    }
    
    // Now do Chi-squared test
    double chisq = 0.0;
    int df = 0;
    for (map<int, int>::iterator d = doubles.begin(); d != doubles.end(); ++d){
        double expected = (double)tot_double*doubfreq[d->first];
        if (expected > 0){
            chisq += pow((double)d->second - expected, 2) / expected;
        }
        ++df;
    }
    // Degrees of freedom = number of categories minus 1
    df -= 1;
    return pchisq(chisq, (double)df);
}

/**
 * Trim the directory off of a full filename path
 */
string filename_nopath(string& filename){
    size_t trim_idx = filename.find_last_of("\\/");
    if (trim_idx != string::npos){
        return filename.substr(trim_idx + 1, filename.length() - trim_idx - 1);
    }
    else{
        return filename;
    }
}


