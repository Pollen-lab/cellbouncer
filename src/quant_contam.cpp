#include <getopt.h>
#include <argp.h>
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
#include <sys/stat.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <zlib.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include <optimML/mixcomp.h>
#include "common.h"
#include "ambient_rna.h"
#include "demux_vcf_io.h"

using std::cout;
using std::endl;
using namespace std;

// ===== Program to profile ambient RNA contamination in cells, =====
//       given output of a demux_vcf run.

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "quant_contam [OPTIONS]\n");
    fprintf(stderr, "Given output of a demux_vcf run, uses the (pre-computed) allele counts\n");
    fprintf(stderr, "to model ambient RNA contamination. Outputs the estimated fraction of\n");
    fprintf(stderr, "each cell's RNA composed of ambient RNA and attempts to find the likeliest\n");
    fprintf(stderr, "mixture of individuals from the VCF that compose the pool of ambient RNA.\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --output_prefix -o The output prefix used in a prior run of demux_vcf\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --num_threads -T The number of parallel threads to use when running\n");
    fprintf(stderr, "       optimization problems (default = 1 = no parallelization)\n");
    fprintf(stderr, "    --bootstrap -b Number of bootstrap replicates to run in order to get\n");
    fprintf(stderr, "       variance on mixture proportions of individuals in the ambient RNA\n");
    fprintf(stderr, "       pool. Default = 100.\n");
    fprintf(stderr, "    --doublet_rate -D Expected probability of doublet droplets (for re-IDing\n");
    fprintf(stderr, "       cells). Default = no expectation. Note that this parameter differs from the\n");
    fprintf(stderr, "       one in demux_vcf: in that program, 0.5 = effectively no prior. In this\n");
    fprintf(stderr, "       program, by default no assumption is made about relative frequencies\n");
    fprintf(stderr, "       of different types of singlets and doublets. If you set this parameter,\n");
    fprintf(stderr, "       however, it will compute the overall frequency of each individual in the\n");
    fprintf(stderr, "       data set (as if bulk), and then use this parameter to determine the expected\n");
    fprintf(stderr, "       frequency of each identity (i.e. if ID1 is 5%% of bulk data, ID2 is 10%%, and\n");
    fprintf(stderr, "       D = 0.1, then expected ID1 singlets are 0.9*0.05 + 0.1*0.05*0.05, ID2 singlets\n");
    fprintf(stderr, "       are 0.9*0.1*0.1 + 0.1*0.1*0.1, and ID1+ID2 doublets are 2*0.1*0.05*0.1.\n");
    fprintf(stderr, "       It will then adjust LLRs to encourage identifying the expected proportion of\n");
    fprintf(stderr, "       each identity. This is most useful in high-contamination data sets where\n");
    fprintf(stderr, "       contamination throws off IDs. If you see many more doublets than expected, \n");
    fprintf(stderr, "       set this parameter; if unsure, ignore.\n");
    fprintf(stderr, "    --run_once -r Standard behavior is to iteratively estimate contam profile\n");
    fprintf(stderr, "       and use it to update cell-individual assignments, then repeat until\n");
    fprintf(stderr, "       log likelihood converges. With this option, it will do this process\n");
    fprintf(stderr, "       once and exit.\n");
    fprintf(stderr, "    --ids -i If you limited the individuals to assign when running demux_vcf\n");
    fprintf(stderr, "       (i.e. your VCF contained extra individuals not in the experiment),\n");
    fprintf(stderr, "       provide the filtered list of individuals here. Should be a text file\n");
    fprintf(stderr, "       with one individual name per line, matching individual names in the VCF.\n");
    fprintf(stderr, "    --ids_doublet -I Similar to --ids/-i argument above, but allows control\n");
    fprintf(stderr, "       over which doublet identities are considered. Here, you can specify\n");
    fprintf(stderr, "       individuals and combinations of two individuals to allow. Doublet\n");
    fprintf(stderr, "       combinations not specified in this file will not be considered.\n");
    fprintf(stderr, "       Single individuals involved in doublet combinations specified here\n");
    fprintf(stderr ,"       but not explicitly listed in the file will still be considered.\n");
    fprintf(stderr, "       Names of individuals must match those in the VCF, and combinations\n");
    fprintf(stderr, "       of two individuals can be specified by giving both individual names\n");
    fprintf(stderr, "       separated by \"+\", with names in either order.\n");
    fprintf(stderr, "    --dump_freqs -d After inferring the ambient RNA profile, write a file containing\n");
    fprintf(stderr, "        alt allele frequencies at each type of site in ambient RNA. This file will\n");
    fprintf(stderr, "        be called [output_prefix].contam.dat\n");
    fprintf(stderr, "    --llr -l Log likelihood ratio cutoff to filter assignments from demux_vcf.\n");
    fprintf(stderr, "        This is the fourth column in the .assignments file. Default = 0 (no filter)\n");
    fprintf(stderr, "    --other_species -s If profiling the ambient RNA is enabled (no -p option),\n");
    fprintf(stderr, "        and your data came from a pool of multiple species, demultiplexed and each\n");
    fprintf(stderr, "        mapped to its species-specific reference genome and then demultiplexed by\n");
    fprintf(stderr, "        individual using within-species SNPs, this option models ambient RNA as a\n");
    fprintf(stderr, "        mixture of all individuals in the VCF, plus RNA from other species.\n");
    fprintf(stderr, "    --error_ref -e The underlying, true rate of misreading reference as\n");
    fprintf(stderr, "        alt alleles (should only reflect sequencing error if variant calls\n");
    fprintf(stderr, "        are reliable; default 0.001)\n");
    fprintf(stderr, "    --error_alt -E The underlying, true rate of misreading alt as reference\n");
    fprintf(stderr, "        alleles (should only reflect sequencing error if variant calls are\n");
    fprintf(stderr, "        reliable; default 0.001)\n"); 
    fprintf(stderr, "    --n_mixprop_trials -N Mixture proportion inference is influenced by initial\n");
    fprintf(stderr, "        guesses. The first time they are inferred, the starting proportions will be\n");
    fprintf(stderr, "        randomly shuffled a number of times equal to this number times the number of\n");
    fprintf(stderr, "        mixture components. Default = 10.\n");
    fprintf(stderr, "    --no_weights -w By default, all observations are weighted by confidence: the log\n");
    fprintf(stderr, "        likelihood ratio of individual ID, divided by the sum of all log likelihood\n");
    fprintf(stderr, "        ratios of assignments of cells to the same individual. This option disables\n");
    fprintf(stderr, "        this weighting and lets all cells contribute equally (although cells with higher\n");
    fprintf(stderr, "        counts will have a stronger influence on the likelihood). You might want to\n");
    fprintf(stderr, "        disable weighting if, for example, you have very unequal numbers of different\n");
    fprintf(stderr, "        individuals and are worried some individual assignments might be mostly noise.\n");
    fprintf(stderr, "        Default behavior would be to give all cells assigned to the noise individual the\n");
    fprintf(stderr, "        same overall weight as all cells assigned to any other individual.\n"); 
    print_libname_help();
    fprintf(stderr, "===== OPTIONAL; FOR INFERRING GENE EXPRESSION =====\n");
    fprintf(stderr, "    --barcodes -B (Optionally gzipped) barcodes file, from MEX-format single cell gene\n");
    fprintf(stderr, "        expression data\n");
    fprintf(stderr, "    --features -F (Optionally gzipped) features file, from MEX-format single cell gene\n");
    fprintf(stderr, "        expression data\n");
    fprintf(stderr, "    --matrix -M (Optionally gizpped) matrix file, from MEX-format single cell gene\n");
    fprintf(stderr, "        expression data\n");
    fprintf(stderr, "    --feature_type -t (OPTIONAL) If --features/-f contains more than one type of data\n");
    fprintf(stderr, "        (i.e. gene expression and feature barcoding), use this to specify which feature\n");
    fprintf(stderr, "        type is RNA-seq (for 10X Genomics, \"Gene Expression\"). By default, includes all\n");
    fprintf(stderr, "        features and does not check.\n");
    fprintf(stderr, "    --clusts -c (RECOMMENDED) cell-cluster assignments computed by another program.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

double ll_gex(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = params[0];
    double p = data_d.at("p");
    return logbinom(n, k, p);
}

void dll_gex(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double n = data_d.at("n");
    double k = params[0];
    double p = data_d.at("p");
    
    // Want derivative wrt k
    results[0] += (digamma(-k+n+1) - digamma(k+1)) + log(p) - log(1.0-p);
}

double ll_gex2(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    double n = data_d.at("n");
    double k = data_d.at("k");
    int i = data_i.at("i");
    double p = data_d.at("p");
    return logbinom(n, k, p*params[i]);
}

void dll_gex2(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
    
    double n = data_d.at("n");
    double k = params[0];
    double p = data_d.at("p");
    int i = data_i.at("i");
    double r = p*params[i];

    double dy_dr = (k-n*r)/(r-r*r);
    results[i] += (dy_dr*p);

}

double ll_multinom(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i){
    
    int num = data_i.at("num");
    double c = data_d.at("c");
    int gi = data_i.at("grp_idx");

    vector<double> p;
    vector<double> n;
    
    char buf[30];
    string bufstr;
    
    double xsum = 1;
    double term2 = 0;
    double term3 = 0;
    double psum = 0.0;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "i_%d", i);
        bufstr = buf;
        int idx = data_i.at(bufstr);
        if (idx == -1){
            break;
        }
        double p_i = c*params[idx] + (1.0-c)*params[gi*num+idx];
        sprintf(&buf[0], "n_%d", i);
        bufstr = buf;
        double n_i = data_d.at(bufstr);
        xsum += n_i;
        term2 += lgammaf(n_i + 1);
        term3 += n_i*log(p_i);
    }
    double term1 = lgammaf(xsum);
    return (term1 - term2 + term3);
}

void dll_multinom(const vector<double>& params,
    const map<string, double>& data_d,
    const map<string, int>& data_i,
    vector<double>& results){
     
    int num = data_i.at("num");
    double c = data_d.at("c");
    double tot = data_d.at("tot");
    int gi = data_i.at("grp_idx");

    char buf[30];
    string bufstr;
    for (int i = 0; i < num; ++i){
        sprintf(&buf[0], "i_%d", i);
        bufstr = buf;
        int idx = data_i.at(bufstr);
        if (idx == -1){
            break;
            //results[idx] -= tot;
            //results[num+idx] -= tot;
        }
        else{
            double p_i = c*params[idx] + (1.0-c)*params[gi*num+idx];
            sprintf(&buf[0], "n_%d", i);
            bufstr = buf;
            double n_i = data_d.at(bufstr);
            double dLL_dpi = (n_i/p_i);
            results[idx] += (dLL_dpi)*c;
            results[gi*num+idx] += (dLL_dpi)*(1.0-c);
        }
    }
}

/**
 * Add in Lagrangian to each element just once at the end
 */
void dll_adjust(const vector<double>& dat_d, const vector<int>& dat_i, vector<double>& G){
    int num = dat_i[0];
    double ncsum = dat_d[0];
    //double sum_totxc = dat_d[0];
    //double sum_totx1mc = dat_d[1];
    for (int i = 0; i < num; ++i){
        //G[i] -= sum_totxc;
        //G[num+i] -= sum_totx1mc;
        G[i] += ncsum;
        G[num+i] -= ncsum;
    }
}

void contam_gex(robin_hood::unordered_map<unsigned long, double>& contam_rate,
    robin_hood::unordered_map<unsigned long, map<int, long int> >& counts, 
    vector<string>& features,
    robin_hood::unordered_map<unsigned long, int>& clusters,
    int n_clusters,
    robin_hood::unordered_map<unsigned long, int>& assn,
    int n_samples,
    map<int, double>& contam_prof,
    int num_threads){
    
    vector<double> grp1;
    vector<vector<double> > grps2;
    vector<vector<double> > ns;
    vector<vector<int> > ns_idx;
    vector<double> tots;
    vector<double> cs;
    vector<int> grp_idx;
    double sum_totxc = 0.0;
    double sum_totx1mc = 0.0;
    double ncsum = 0.0;
    
    vector<double> grps2tot;
    int n_grp2 = n_clusters;
    if (n_clusters == 0){
        n_grp2 = 1;
    }
    for (int i = 0; i < n_grp2; ++i){
        vector<double> v;
        grps2.push_back(v);
        grps2tot.push_back((double)features.size());
    }

    for (int i = 0; i < features.size(); ++i){
        grp1.push_back(1.0);
        //grp1.push_back(1.0/(double)features.size());
        for (int j = 0; j < n_grp2; ++j){
            grps2[j].push_back(1.0);
        }
        vector<double> v;
        ns.push_back(v);
        vector<int> v2;
        ns_idx.push_back(v2);
    }
    double grp1tot = (double)features.size();
    int num_added = 0;
    double grp2tot = (double)features.size();

    for (robin_hood::unordered_map<unsigned long, double>::iterator cr = contam_rate.begin();
        cr != contam_rate.end(); ++cr){

        int clust = 0;
        if (clusters.size() > 0){
            if (clusters.count(cr->first) == 0){
                clust = -1;
            }
            else{
                clust = clusters[cr->first];
            }
        }
        int cell_assn = -1;
        if (assn.count(cr->first) > 0){
            cell_assn = assn[cr->first];
        }
        if (clust < 0 && cell_assn < -1){
            continue;
        }

        double celltot = 0;
        
        int i = 0;
        int cprev = -1;
        for (map<int, long int>::iterator c = counts[cr->first].begin(); c != counts[cr->first].end(); ++c){
            double count = (double)c->second;
            
            if (cell_assn >= 0){
                if (cell_assn >= n_samples){
                    pair<int, int> comb = idx_to_hap_comb(cell_assn, n_samples);
                    grp1[comb.first] += 0.5*contam_prof[comb.first]*count;
                    grp1[comb.second] += 0.5*contam_prof[comb.second]*count;
                    grp1tot += 0.5*contam_prof[comb.first]*count + 0.5*contam_prof[comb.second]*count;
                }
                else{ 
                    grp1[c->first] += contam_prof[cell_assn] * count;
                    grp1tot += contam_prof[cell_assn] * count;         
                }
            }
            if (clust >= 0){
                grps2[clust][c->first] += count;
                grps2tot[clust] += count;

                celltot += count;
                ns_idx[i].push_back(c->first);
                for (int j = cprev+1; j < c->first; ++j){
                    ns[j].push_back(0);
                }
                ns[c->first].push_back(count);
                i++;
                cprev = c->first;
            }
        }
        if (clust < 0){
            continue;
        }
        for (int j = i; j < features.size(); ++j){
            ns_idx[j].push_back(-1);
        }
        for (int j = cprev + 1; j < features.size(); ++j){
            ns[j].push_back(0.0);
        }
        tots.push_back(celltot);
        cs.push_back(cr->second);
        grp_idx.push_back(clust + 1);
        //sum_totxc += (celltot*cr->second);
        //sum_totx1mc += (celltot*(1.0-cr->second));
        //ncsum += celltot*cr->second;

        ++num_added;
    }
    
    for (int j = 0; j < grps2.size(); ++j){ 
        for (int i = 0; i < features.size(); ++i){
            grps2[j][i] /= grp2tot;
            if (j == 0){
                grp1[i] /= grp1tot;
            }
        }
    }
    vector<double> paramsx = {};
    optimML::multivar_ml_solver mnsolver(paramsx, ll_multinom, dll_multinom);
    if (num_threads > 1){
        mnsolver.set_threads(num_threads);
        mnsolver.set_bfgs_threads(num_threads);
    }
    mnsolver.add_param_grp(grp1);
    for (int i = 0; i < grps2.size(); ++i){
        mnsolver.add_param_grp(grps2[i]);
    }
    mnsolver.add_data_fixed("num", (int)features.size());
    mnsolver.add_data("tot", tots);
    mnsolver.add_data("c", cs);
    mnsolver.add_data("grp_idx", grp_idx);
    char buf[30];
    for (int i = 0; i < features.size(); ++i){
        sprintf(&buf[0], "n_%d", i);
        string bufstr = buf;
        mnsolver.add_data(bufstr, ns[i]);
        sprintf(&buf[0], "i_%d", i);
        bufstr = buf;
        mnsolver.add_data(bufstr, ns_idx[i]);
    }
    //vector<double> hookparamsd{ ncsum };
    //vector<int> hookparamsi{ (int)features.size() };
    //mnsolver.add_gradient_hook(dll_adjust, hookparamsd, hookparamsi);
    //mnsolver.set_maxiter(10000);
    mnsolver.set_delta(1);
    fprintf(stderr, "Inferring ambient RNA expression profile...\n");
    mnsolver.solve();
    double sum1 = 0.0;
    vector<double> sums2;
    for (int i = 0; i < n_grp2; ++i){
        sums2.push_back(0.0);
    }
    for (int i = 0; i < features.size(); ++i){
        fprintf(stdout, "%s\t%e", features[i].c_str(), mnsolver.results[i]);
        sum1 += mnsolver.results[i];
        for (int j = 0; j < n_grp2; ++j){
            fprintf(stdout, "\t%e", mnsolver.results[features.size()*(j+1) + i]);
            sums2[j] += mnsolver.results[features.size()*(j+1) + i];
        }
        fprintf(stdout, "\n");
    }
    fprintf(stderr, "sum1 %f\n", sum1);
    fprintf(stderr, "sums2:\n");
    for (int j = 0; j < n_grp2; ++j){
        fprintf(stderr, "  %f\n", sums2[j]);
    }
    exit(0);

    // Model: Binomial(cell_tot, sum(contam_vec_gene*count_gene), contam_rate)
    // (n,k,p)

    //vector<vector<double> > mixfracs;
    //vector<map<int, double> > mixfracs;
    //vector<double> n;
    //vector<double> k;
    //vector<double> p;
    
    /* 
    vector<int> gene_idx;
    for (robin_hood::unordered_map<unsigned long, double>::iterator cr = contam_rate.begin();
        cr != contam_rate.end(); ++cr){
        long int celltot = 0;
        
        for (map<int, long int>::iterator c = counts[cr->first].begin(); c != 
            counts[cr->first].end(); ++c){
            gene_idx.push_back(c->first);
            k.push_back((double)c->second);    
            p.push_back(cr->second);
            celltot += c->second;
        }
        for (int i = 0; i < counts[cr->first].size(); ++i){
            n.push_back((double)celltot);
        }
    }
    vector<double> params1 = {};
    optimML::multivar_ml_solver solver1(params1, ll_gex2, dll_gex2);
    vector<double> grp;
    for (int i = 0; i < features.size(); ++i){
        grp.push_back(1.0/(double)features.size());
    } 
    solver1.add_param_grp(grp);
    solver1.add_data("n", n);
    solver1.add_data("k", k);
    solver1.add_data("p", p);
    solver1.add_data("i", gene_idx);
    fprintf(stderr, "solving\n");
    solver1.solve();
    fprintf(stderr, "LL %f\n", solver1.log_likelihood);
    for (int i = 0; i < solver1.results.size(); ++i){
        fprintf(stdout, "%s\t%f\n", features[i].c_str(), solver1.results[i]);
    }
    exit(0);
    */

    /*
    for (robin_hood::unordered_map<unsigned long, double>::iterator cr = contam_rate.begin();
        cr != contam_rate.end(); ++cr){
        long int celltot = 0;
        map<int, double> row;    
        for (map<int, long int>::iterator c = counts[cr->first].begin(); c !=
            counts[cr->first].end(); ++c){
            celltot += c->second;
            row.insert(make_pair(c->first, (double)c->second));
        }
        mixfracs.push_back(row);
        n.push_back((double)celltot);
        p.push_back(cr->second);
    }
    
    fprintf(stderr, "gathered %ld %ld %ld\n", n.size(), p.size(), mixfracs.size());

    vector<double> params = { };
    optimML::multivar_ml_solver solver(params, ll_gex, dll_gex);
    solver.add_mixcomp(mixfracs, features.size());
    solver.add_data("n", n);
    solver.add_data("p", p);
    fprintf(stderr, "solving\n");
    solver.solve();
    fprintf(stderr, "LL %f\n", solver.log_likelihood);
    
    vector<double> genesum_before;
    vector<double> genesum_after;
    for (int i = 0; i < solver.results_mixcomp.size(); ++i){
        genesum_before.push_back(0.0);
        genesum_after.push_back(0.0);
        fprintf(stdout, "%s\t%f\n", features[i].c_str(), solver.results_mixcomp[i]);
    }
    exit(0);
    for (robin_hood::unordered_map<unsigned long, double>::iterator c = contam_rate.begin(); c != 
        contam_rate.end(); ++c){
        for (map<int, long int>::iterator count = counts[c->first].begin(); count != counts[c->first].end();
            ++count){
            genesum_before[count->first] += (double)count->second;
            double count_adj = (double)count->second;
            double contam_prop = c->second * solver.results_mixcomp[count->first];
            count_adj *= (1.0 - contam_prop);
            genesum_after[count->first] += count_adj;
        }
    }

    for (int i = 0; i < solver.results_mixcomp.size(); ++i){
        //fprintf(stdout, "%s\t%f\n", features[i].c_str(), solver.results_mixcomp[i]);
        fprintf(stdout, "%s\t%f\t%f\n", features[i].c_str(), genesum_before[i], genesum_after[i]); 
    }
    */
}

int parse_clustfile(string& filename, 
    robin_hood::unordered_map<unsigned long, int>& clusts){
    
    ifstream inf(filename.c_str());
    string bcstr;
    string clustn;
    set<string> clustnames;
    robin_hood::unordered_map<unsigned long, string> ctmp;
    while (inf >> bcstr >> clustn){
        clustnames.insert(clustn);
        unsigned long ul = bc_ul(bcstr);
        ctmp.emplace(ul, clustn);
    }

    map<string, int> clust2idx;
    int idx = 0;
    for (set<string>::iterator s = clustnames.begin(); s != clustnames.end(); ++s){
        clust2idx.insert(make_pair(*s, idx));
        idx++;
    }
    for (robin_hood::unordered_map<unsigned long, string>::iterator t = ctmp.begin(); 
        t != ctmp.end(); ++t){
        clusts.emplace(t->first, clust2idx[t->second]);
    }
    return clustnames.size();
}

int main(int argc, char *argv[]) {    

    static struct option long_options[] = {
       {"output_prefix", required_argument, 0, 'o'},
       {"other_species", required_argument, 0, 's'},
       {"error_ref", required_argument, 0, 'e'},
       {"error_alt", required_argument, 0, 'E'},
       {"doublet_rate", required_argument, 0, 'D'},
       {"llr", required_argument, 0, 'l'},
       {"n_mixprop_trials", required_argument, 0, 'N'},
       {"no_weights", no_argument, 0, 'w'},
       {"dump_freqs", no_argument, 0, 'd'},
       {"ids", required_argument, 0, 'i'},
       {"ids_doublet", required_argument, 0, 'I'},
       {"libname", required_argument, 0, 'n'},
       {"cellranger", no_argument, 0, 'C'},
       {"seurat", no_argument, 0, 'S'},
       {"underscore", no_argument, 0, 'U'},
       {"run_once", no_argument, 0, 'r'},
       {"bootstrap", required_argument, 0, 'b'},
       {"barcodes", required_argument, 0, 'B'},
       {"features", required_argument, 0, 'F'},
       {"matrix", required_argument, 0, 'M'},
       {"feature_type", required_argument, 0, 't'},
       {"clusts", required_argument, 0, 'c'},
       {"num_threads", required_argument, 0, 'T'},
       {0, 0, 0, 0} 
       
    };
    
    // Set default values
    string output_prefix = "";
    bool inter_species = false;
    double error_ref = 0.001;
    double error_alt = 0.001;
    double llr = 0.0;
    int n_mixprop_trials = 10;
    bool weight = true;
    bool dump_freqs = false;
    string idfile;
    bool idfile_given = false;
    string idfile_doublet;
    bool idfile_doublet_given = false;
    string libname = "";
    bool cellranger = false;
    bool seurat = false;
    bool underscore = false;
    bool run_once = false;
    int bootstrap = 100;
    double doublet_rate = -1.0;
    int num_threads = 0;

    string barcodesfile = "";
    string featuresfile = "";
    string matrixfile = "";
    string feature_type = "";
    string clustfile = "";

    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
    while((ch = getopt_long(argc, argv, "o:e:E:l:N:i:I:n:b:D:B:F:M:t:c:T:rsCSUdwh", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'i':
                idfile = optarg;
                idfile_given = true;
                break;
            case 'I':
                idfile_doublet = optarg;
                idfile_doublet_given = true;
                break;
            case 's':
                inter_species = true;
                break;
            case 'D':
                doublet_rate = atof(optarg);
                break;
            case 'e':
                error_ref = atof(optarg);
                break;
            case 'E':
                error_alt = atof(optarg);
                break;
            case 'l':
                llr = atof(optarg);
                break;
            case 'N':
                n_mixprop_trials = atoi(optarg);
                break;
            case 'w':
                weight = false;
                break;
            case 'd':
                dump_freqs = true;
                break;
            case 'n':
                libname = optarg;
                break;
            case 'C':
                cellranger = true;
                break;
            case 'S':
                seurat = true;
                break;
            case 'U':
                underscore = true;
                break;
            case 'r':
                run_once = true;
                break;
            case 'b':
                bootstrap = atoi(optarg);
                break;
            case 'B':
                barcodesfile = optarg;
                break;
            case 'F':
                featuresfile = optarg;
                break;
            case 'M':
                matrixfile = optarg;
                break;
            case 't':
                feature_type = optarg;
                break;
            case 'c':
                clustfile = optarg;
                break;
            case 'T':
                num_threads = atoi(optarg);
                break;
            default:
                help(0);
                break;
        }    
    }
     
    // Error check arguments.
    if (output_prefix == ""){
        fprintf(stderr, "ERROR: output_prefix required\n");
        exit(1);
    }
    if (error_ref <= 0 || error_ref >= 1.0 || error_alt <= 0 || error_alt >= 1.0){
        fprintf(stderr, "ERROR: error rates must be between 0 and 1, exclusive.\n");
        exit(1);
    }
    if (n_mixprop_trials < 0){
        fprintf(stderr, "ERROR: --n_mixprop_trials must be >= 0\n");
        exit(1);
    }
    if (idfile_given && idfile_doublet_given){
        fprintf(stderr, "ERROR: only one of -i and -I is allowed.\n");
        exit(1);
    }
    if (bootstrap <= 0){
        fprintf(stderr, "WARNING: bootstrapping disabled. Ambient RNA pool proportions will \
be reported without concentration parameters (variance will be unknown).\n");
    }
    if (doublet_rate != -1 && (doublet_rate < 0 || doublet_rate > 1)){
        fprintf(stderr, "ERROR: --doublet_rate/-D must be between 0 and 1, inclusive.\n");
        exit(1);
    }
    if ((barcodesfile != "" || featuresfile != "" || matrixfile != "") &&
        !(barcodesfile != "" && featuresfile != "" && matrixfile != "")){
        fprintf(stderr, "ERROR: if inferring gene expression profile, you must provide all\n");
        fprintf(stderr, "three of --barcodes/-B, --features/-F, and --matrix/-M\n");
        exit(1);
    }
    if (barcodesfile != "" && clustfile == ""){
        fprintf(stderr, "WARNING: inferring expression profile of contamination without \
cluster information. Assuming one default expression profile for each individual (results will \
be inaccurate if there is much cell type heterogeneity).\n");
    }
    if (clustfile != "" && barcodesfile == ""){
        fprintf(stderr, "ERROR: --clusters/-c only applicable when loading gene expression data\n");
        exit(1);
    }

    // Load files
    string sample_name = output_prefix + ".samples";
    vector<string> samples;
    if (file_exists(sample_name)){
        load_samples(sample_name, samples);
    }
    else{
        fprintf(stderr, "ERROR: no samples file found for %s. Please run demux_vcf with\n", 
            output_prefix.c_str());
        fprintf(stderr, "same output prefix.\n");
        exit(1);
    }
    
    map<pair<int, int>, map<int, float> > exp_match_fracs;
    string expfrac_name = output_prefix + ".condf";
    if (file_exists(expfrac_name)){
        load_exp_fracs(expfrac_name, exp_match_fracs);
    }
    else{
        fprintf(stderr, "ERROR: no conditional matching probability file found for %s.\n", 
            output_prefix.c_str());
        fprintf(stderr, "Please re-run demux_vcf with the same VCF file and output prefix, \
but specify the -F option to create this file. Then re-run this program.\n");
        exit(1); 
    }

    set<int> allowed_ids;
    set<int> allowed_ids2;

    if (idfile_given){
        parse_idfile(idfile, samples, allowed_ids, allowed_ids2, true);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile.c_str());
        }
    }
    if (idfile_doublet_given){
        parse_idfile(idfile_doublet, samples, allowed_ids, allowed_ids2, false);
        if (allowed_ids.size() == 0){
            fprintf(stderr, "No valid individual names found in file %s; allowing \
all possible individuals\n", idfile_doublet.c_str());
        }
    }

    // 1 thread means 0 threads (don't launch extra processes)
    if (num_threads <= 1){
        num_threads = 0;
    }
    
    // Map cell barcodes to numeric IDs of best individual assignments 
    robin_hood::unordered_map<unsigned long, int> assn;

    // Map cell barcodes to log likelihood ratio of best individual assignments
    robin_hood::unordered_map<unsigned long, double> assn_llr;
    
    string assn_name = output_prefix + ".assignments";
    if (file_exists(assn_name)){
        fprintf(stderr, "Loading assignments...\n");
        load_assignments_from_file(assn_name, assn, assn_llr, samples);
        if (llr > 0.0){
            // Filter assignments.
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end();){
                if (assn_llr[a->first] <= llr){
                    assn_llr.erase(a->first);
                    // Can't erase while iterating the normal way
                    // https://github.com/martinus/robin-hood-hashing/issues/18
                    a = assn.erase(a);
                }
                else{
                    ++a;
                }
            }
            if (assn.size() == 0){
                fprintf(stderr, "ERROR: LLR filter too high; no assignments left to use.\n");
                exit(1);
            }
        }
    }
    else{
        fprintf(stderr, "ERROR: no assignments found for %s. Please run demux_vcf with same\n", 
            output_prefix.c_str());
        fprintf(stderr, "output prefix.\n");
        exit(1); 
    }
    
    // Store allele counts
    robin_hood::unordered_map<unsigned long, map<pair<int, int>,
        map<pair<int, int>, pair<float, float> > > > indv_allelecounts;
    string counts_name = output_prefix + ".counts";
    if (file_exists(counts_name)){
        fprintf(stderr, "Loading counts...\n");
        load_counts_from_file(indv_allelecounts, samples, counts_name, allowed_ids); 
    }
    else{
        fprintf(stderr, "ERROR: no counts found for %s. Please run demux_vcf with same\n",
            output_prefix.c_str());
        fprintf(stderr, "output prefix.\n");
        exit(1);
    }
    
    double llprev = 0.0;
    double delta = 999;
    double delta_thresh = 0.1;
    
    map<int, double> contam_prof;
    map<int, double> contam_prof_conc;
    robin_hood::unordered_map<unsigned long, double> contam_rate;
    robin_hood::unordered_map<unsigned long, double> contam_rate_se;
    int nits = 0;
    while (delta > delta_thresh){
        fprintf(stderr, "===== ITERATION %d =====\n", nits+1);
        contamFinder cf(indv_allelecounts, assn, assn_llr, exp_match_fracs, samples.size(),
            allowed_ids, allowed_ids2);
        cf.set_doublet_rate(doublet_rate);
        cf.set_num_threads(num_threads);
        // Initialize to whatever was the final estimate last time 
        if (nits > 0){
            cf.set_init_contam_prof(contam_prof);
        }
        if (nits > 0){
            double meanc = 0.0;
            double cfrac = 1.0/(double)contam_rate.size();
            for (robin_hood::unordered_map<unsigned long, double>::iterator c = contam_rate.begin();
                c != contam_rate.end(); ++c){
                meanc += cfrac * c->second;
            }
            cf.set_init_c(meanc);
        } 
        // Do standard initialization
        cf.set_error_rates(error_ref, error_alt);
        if (inter_species){
            cf.model_other_species();
        } 
        cf.set_mixprop_trials(n_mixprop_trials);
        if (weight){
            cf.use_weights();
        }
        cf.fit(); 
        double ll = cf.compute_ll();
        
        if (run_once){
            // Break out of cycle
            assn = cf.assn;
            assn_llr = cf.assn_llr;
            contam_prof = cf.contam_prof;
            contam_rate = cf.contam_rate;
            contam_rate_se = cf.contam_rate_se;
            
            delta = 0;
        }
        else{
            if (llprev == 0 || ll > llprev){
                assn = cf.assn;
                assn_llr = cf.assn_llr;
                contam_prof = cf.contam_prof;
                contam_rate = cf.contam_rate;
                contam_rate_se = cf.contam_rate_se;
            }
            fprintf(stderr, " -- Log likelihood: %f", ll);
            if (llprev != 0){
                delta = ll-llprev;
                fprintf(stderr, " delta = %f\n", delta);
            }
            else{
                fprintf(stderr, "\n");
            }
            llprev = ll;
            nits++;
        }
        if (delta <= delta_thresh && bootstrap > 0){
            // Do bootstrapping
            cf.assn = assn;
            cf.assn_llr = assn_llr;
            cf.contam_prof = contam_prof;
            cf.contam_rate = contam_rate;
            cf.contam_rate_se = contam_rate_se;
            fprintf(stderr, "Computing Dirichlet concentration parameters \
on mixture proportions...\n");
            cf.bootstrap_amb_prof(bootstrap, contam_prof_conc);
        }
    }

    // Write contamination profile to disk
    {
        string fname = output_prefix + ".contam_prof";
        FILE* outf = fopen(fname.c_str(), "w");
        fprintf(stderr, "Writing contamination profile to disk...\n");
        dump_contam_prof(outf, contam_prof, contam_prof_conc, samples);
        fclose(outf);
    }
    // Write contamination rate (and standard error) per cell to disk
    {
        string fname = output_prefix + ".contam_rate";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_contam_rates(outf, contam_rate, contam_rate_se, samples,
            libname, cellranger, seurat, underscore);
        fclose(outf);
    }
    /*
    if (dump_freqs){
        // Write alt allele matching frequencies to a file
        string fname = output_prefix + ".contam.dat";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_amb_fracs(outf, cf.amb_mu);
        fclose(outf);
    }
    */
    // Write refined assignments to disk
    {
        string fname = output_prefix + ".decontam.assignments";
        FILE* outf = fopen(fname.c_str(), "w");
        dump_assignments(outf, assn, assn_llr, samples, libname, 
            cellranger, seurat, underscore);
        fclose(outf); 
    }
    if (barcodesfile != "" && featuresfile != "" && matrixfile != ""){
        robin_hood::unordered_map<unsigned long, map<int, long int> > mtx;
        vector<string> features;
        fprintf(stderr, "Loading gene expression data...\n");
        parse_mex(barcodesfile, featuresfile, matrixfile, mtx, features, feature_type);
        robin_hood::unordered_map<unsigned long, int> clusts;
        int nclusts = 0;
        if (clustfile != ""){
            nclusts = parse_clustfile(clustfile, clusts);
        }
        else{
            fprintf(stderr, "Using cell identities as clusters\n");
            for (robin_hood::unordered_map<unsigned long, int>::iterator a = assn.begin();
                a != assn.end(); ++a){
                if (a->second < samples.size()){
                    clusts.emplace(a->first, a->second);
                }
            } 
            nclusts = samples.size();
        }
        contam_gex(contam_rate, mtx, features, clusts, nclusts, assn, samples.size(), contam_prof, num_threads); 
    }
    return 0;
}
