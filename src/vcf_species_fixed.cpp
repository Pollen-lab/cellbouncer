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
#include <random>
#include <functional>
#include <utility>
#include <math.h>
#include <float.h>
#include <htslib/vcf.h>
#include <htswrapper/gzreader.h>
#include <zlib.h>

using std::cout;
using std::endl;
using namespace std;

/**
 * Print a help message to the terminal and exit.
 */
void help(int code){
    fprintf(stderr, "vcf_species_fixed [OPTIONS]\n");
    fprintf(stderr, "Given a VCF file containing variatns called across multiple species,\n");
    fprintf(stderr, "along with a file mapping individuals to species, finds sites at which\n");
    fprintf(stderr, "all species are fixed for different alleles, and outputs a new VCF\n");
    fprintf(stderr, "with one entry per species, showing only fixed sites.\n");
    fprintf(stderr, "Additionally, you can include F1s in the output file as follows:\n");
    fprintf(stderr, " --f1 Chinobo,Chimp,Bonobo\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v A VCF/BCF file listing variants.\n");
    fprintf(stderr, "    --species -s A file mapping individual ID to species, one per line,\n");
    fprintf(stderr, "      tab-separated.\n");
    fprintf(stderr, "    --f1 -f Simulate F1 hybrid between species in output. You can exclude\n");
    fprintf(stderr, "      this option or supply it multiple times. Syntax is name,species1,species2\n");
    fprintf(stderr, "      In other words, to include a column in the output VCF called Hybrid, which\n");
    fprintf(stderr, "      will have the expected genotypes of a hybrid between species1 and species2\n");
    fprintf(stderr, "      you can supply --f1 Hybrid,species1,species2.\n");
    fprintf(stderr, "    --output -o Output file name.\n");
    fprintf(stderr, "    --output_fmt -O Output file format (from HTSLib: v = VCF; z = gzVCF;\n");
    fprintf(stderr, "      b = BCF; etc. Default: z\n");
    fprintf(stderr, "    --help -h Display this message and exit.\n");
    exit(code);
}

/**
 * Returns true if fixed between species
 * false if not fixed between species
 */
bool proc_bcf_record(bcf1_t* bcf_record,
    bcf_hdr_t* bcf_header,
    int num_samples,
    int num_species,
    vector<int>& sample_species,
    vector<int>& species_ref,
    vector<int>& species_alt){
    
    species_ref.clear();
    species_alt.clear();
    for (int i = 0; i < num_species; ++i){
        species_ref.push_back(0);
        species_alt.push_back(0);
    }

    // Load ref/alt alleles and other stuff
    // This puts alleles in bcf_record->d.allele[index]
    // Options for parameter 2:
    
    // BCF_UN_STR  1       // up to ALT inclusive
    // BCF_UN_FLT  2       // up to FILTER
    // BCF_UN_INFO 4       // up to INFO
    // BCF_UN_SHR  (BCF_UN_STR|BCF_UN_FLT|BCF_UN_INFO) // all shared information
    // BCF_UN_FMT  8                           // unpack format and each sample
    // BCF_UN_IND  BCF_UN_FMT                  // a synonymo of BCF_UN_FMT
    // BCF_UN_ALL (BCF_UN_SHR|BCF_UN_FMT) // everything
    
    bcf_unpack(bcf_record, BCF_UN_STR);
    bool pass = true;
    // pass = biallelic, no indels, ref/alt both A/C/G/T
    for (int i = 0; i < bcf_record->n_allele; ++i){
        if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
            strcmp(bcf_record->d.allele[i], "C") != 0 &&
            strcmp(bcf_record->d.allele[i], "G") != 0 && 
            strcmp(bcf_record->d.allele[i], "T") != 0){
            pass = false;
            break;
        }
    }
    if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
        pass = false;
    }
    if (pass){
        // Get all available genotypes.
        int32_t* gts = NULL;
        int n_gts = 0;
        int nmiss = 0;
        int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
        if (num_loaded <= 0){
            fprintf(stderr, "ERROR loading genotypes at %s %ld\n", 
                bcf_hdr_id2name(bcf_header, bcf_record->rid), (long int) bcf_record->pos);
            exit(1);
        }
        
        // Assume ploidy = 2
        int ploidy = 2;
        //int ploidy = n_gts / num_samples; 
       
        int nref = 0;
        int nalt = 0;
         
        for (int i = 0; i < num_samples; ++i){
            int32_t* gtptr = gts + i*ploidy;
            int species_idx = sample_species[i];
            if (species_idx >= 0){   
                if (!bcf_gt_is_missing(gtptr[0])){
                    if (bcf_gt_allele(gtptr[0]) == 0){
                        nref++;        
                        species_ref[species_idx]++;
                    }
                    else{
                        nalt++;
                        species_alt[species_idx]++;
                    }
                    if (bcf_gt_allele(gtptr[1]) == 0){
                        nref++;
                        species_ref[species_idx]++;
                    }
                    else{
                        nalt++;
                        species_alt[species_idx]++;
                    }
                }
                if (species_ref[species_idx] > 0 && species_alt[species_idx] > 0){
                    return false;
                }
            }
        }
        // Exclude sites without polymorphism
        if (nref == 0 || nalt == 0){
            return false;
        } 
        // At this point, there is no within-species polymorphism, and there is polymorphism.
        // Therefore, this must be a fixed difference across species.
        return true;
    }
    // Site did not pass filters.
    return false;
}

void parse_species(string speciesfile, map<string, string>& id2species){
    gzreader reader(speciesfile);
    reader.delimited('\t');
    /*
    inf = ifstream(speciesfile);
    string id;
    string species;
    while (inf >> id >> species){
    */
    while(reader.next()){
        string id = reader.fields[0];
        string species = reader.fields[1];
        if (id2species.count(id) > 0 && id2species[id] != species){
            fprintf(stderr, "ERROR: id %s mapped to multiple species in %s\n", 
                id.c_str(), speciesfile.c_str());
            exit(1);
        }
        id2species.insert(make_pair(id, species));
    }
}

void parse_f1str(string f1str, vector<string>& names, vector<string>& s1, vector<string>& s2){
    string part1;
    string part2;
    string part3;
    size_t pos = f1str.find(",");
    if (pos != string::npos){
        part1 = f1str.substr(0, pos);
        f1str = f1str.substr(pos+1, f1str.length()-pos-1);
        pos = f1str.find(",");
        if (pos != string::npos){
            part2 = f1str.substr(0, pos);
            part3 = f1str.substr(pos+1, f1str.length()-pos-1);
            names.push_back(part1);
            s1.push_back(part2);
            s2.push_back(part3);
            return;
        }
    }
    fprintf(stderr, "ERROR parsing F1 identity string %s\n", f1str.c_str());
    fprintf(stderr, "Expected format: name,species1,species2\n");
    exit(1);
}

int main(int argc, char *argv[]) {    
    
    static struct option long_options[] = {
       {"vcf", required_argument, 0, 'v'},
       {"species", required_argument, 0, 's'},
       {"f1", required_argument, 0, 'f'},
       {"output", required_argument, 0, 'o'},
       {"output_fmt", required_argument, 0, 'O'},
       {0, 0, 0, 0} 
    };
    
    // Set default values
    string vcf_file = "";
    string outfile = "";
    string speciesfile = "";
    string outfmt = "z";
    vector<string> f1_names;
    vector<string> f1_s1;
    vector<string> f1_s2;
    int option_index = 0;
    int ch;
    
    if (argc == 1){
        help(0);
    }
   while((ch = getopt_long(argc, argv, "v:o:s:O:f:h", long_options, &option_index )) != -1){
        switch(ch){
            case 0:
                // This option set a flag. No need to do anything here.
                break;
            case 'h':
                help(0);
                break;
            case 'v':
                vcf_file = optarg;
                break;
            case 'f':
                parse_f1str(optarg, f1_names, f1_s1, f1_s2);
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'O':
                outfmt = optarg;
                break;
            case 's':
                speciesfile = optarg;
                break;
            default:
                help(0);
                break;
        }    
    }
    
    // Error check arguments.
    if (vcf_file == ""){
        fprintf(stderr, "ERROR: VCF file required\n");
        exit(1);
    }
    if (speciesfile == ""){
        fprintf(stderr, "ERROR: species file required\n");
        exit(1);
    }
    if (outfile == ""){
        fprintf(stderr, "ERROR: output file name required\n");
        exit(1);
    }

    if (outfmt == "z"){
        // Clean up output file name
        if (outfile.rfind(".vcf") == string::npos){
            outfile += ".vcf.gz";
        }
        else if (outfile.rfind(".gz") == string::npos){
            outfile += ".gz";
        }
    }
    
    // Read individual -> species assignments
    map<string, string> id2species;
    parse_species(speciesfile, id2species);
    
    // Store the names of all individuals in the VCF
    vector<string> samples;
    vector<string> samples_species;

    // First pass: count occurrences of each branch, only considering sites without
    // missing genotypes.
    bcf_hdr_t* bcf_header;
    bcf1_t* bcf_record = bcf_init();
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR interpreting %s as BCF format.\n", vcf_file.c_str());
        exit(1);
    }
    bcf_header = bcf_hdr_read(bcf_reader);
    int num_samples = bcf_hdr_nsamples(bcf_header);
    for (int i = 0; i < num_samples; ++i){
        samples.push_back(bcf_header->samples[i]);
        if (id2species.count(bcf_header->samples[i]) == 0){
            fprintf(stderr, "WARNING: sample %s not assigned to a species\n", bcf_header->samples[i]);
        }
    }
    
    // Write gz-compressed by default
    string writestr = "w" + outfmt;
    htsFile* outf = hts_open(outfile.c_str(), writestr.c_str());
    

    // Convert sample names to species indices
    vector<int> sample_species;
    map<string, int> species2idx;
    vector<string> idx2species;
    int species_idx = 0;
    for (vector<string>::iterator samp = samples.begin(); samp != samples.end(); ++samp){
        if (id2species.count(*samp) == 0){
            sample_species.push_back(-1);
        }
        else{
            string spec = id2species[*samp];
            int idx;
            if (species2idx.count(spec) > 0){
                idx = species2idx[spec];
            }
            else{
                idx = species_idx;
                species2idx.insert(make_pair(spec, idx));
                species_idx++;
                idx2species.push_back(spec);
            }
            sample_species.push_back(idx);
        }
    }
    
    // Create an output copy of header
    bcf_hdr_t *out_hdr = bcf_hdr_dup(bcf_header);
    
    // Drop old sample names from output header
    int ret = bcf_hdr_set_samples(out_hdr, NULL, 0);
    
    // Add new sample names (species)
    for (int i = 0; i < idx2species.size(); ++i){
        bcf_hdr_add_sample(out_hdr, idx2species[i].c_str());
    }

    int num_species = species_idx;
    vector<pair<int, int> > f1_species_idx;
    for (int i = 0; i < f1_names.size(); ++i){
        int i1;
        int i2;
        if (species2idx.count(f1_s1[i]) == 0){
            fprintf(stderr, "ERROR: unknown species name in F1 combination: %s\n", f1_s1[i].c_str());
            exit(1);
        }    
        else{
            i1 = species2idx[f1_s1[i]];
        }
        if (species2idx.count(f1_s2[i]) == 0){
            fprintf(stderr, "ERROR: unknown species name in F1 combination: %s\n", f1_s2[i].c_str());
            exit(1);
        }
        else{
            i2 = species2idx[f1_s2[i]];
        }
        bcf_hdr_add_sample(out_hdr, f1_names[i].c_str());
        f1_species_idx.push_back(make_pair(i1, i2));
    } 

    ret = bcf_hdr_sync(out_hdr);

    // Write out new header
    int write_success = bcf_hdr_write(outf, out_hdr);

    long int kept = 0;
    
    long int nsnp = 0;

    vector<int> species_ref;
    vector<int> species_alt;
    
    vector<int32_t> species_gt((num_species + f1_s1.size()) * 2);
    
    bcf1_t* out_rec = bcf_init();
    

    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        
        if (bcf_record->n_allele == 2){ 
            bool pass;        
            
            if (proc_bcf_record(bcf_record, bcf_header, num_samples,
                num_species, sample_species, species_ref, species_alt)){
                
                // This is a fixed difference.
                // Copy the record and set it up to point to species-fixed genotypes.
                bcf_copy(out_rec, bcf_record);
                bcf_unpack(out_rec, BCF_UN_ALL);
                // Remove all FORMAT fields
                for (int i = out_rec->n_fmt-1; i >= 0; i--) {
                    bcf_fmt_t *fmt = &out_rec->d.fmt[i];
                    const char *key = bcf_hdr_int2id(bcf_header, BCF_DT_ID, fmt->id);
                    bcf_update_format(bcf_header, out_rec, key, NULL, 0, BCF_HT_INT);
                }
                // Add a new genotype per species (and simulated F1s)
                for (int i = 0; i < species_ref.size(); ++i){
                    if (species_ref[i] > 0){
                        // Ref allele
                        species_gt[i * 2] = bcf_gt_unphased(0);
                        species_gt[i * 2 + 1] = bcf_gt_unphased(0);
                    }
                    else{
                        // Alt allele
                        species_gt[i * 2] = bcf_gt_unphased(1);
                        species_gt[i * 2 + 1] = bcf_gt_unphased(1);
                    }
                }
                for (int i = 0; i < f1_species_idx.size(); ++i){
                    int s1i = f1_species_idx[i].first;
                    int s2i = f1_species_idx[i].second;
                    int j = num_species + i;
                    if (species_ref[s1i] > 0){
                        if (species_ref[s2i] > 0){
                            species_gt[j * 2] = bcf_gt_unphased(0);
                            species_gt[j * 2 + 1] = bcf_gt_unphased(0);
                        }
                        else{
                            species_gt[j * 2] = bcf_gt_unphased(0);
                            species_gt[j * 2 + 1] = bcf_gt_unphased(1);
                        }
                    }
                    else{
                        if (species_ref[s2i] > 0){
                            species_gt[j * 2] = bcf_gt_unphased(0);
                            species_gt[j * 2 + 1] = bcf_gt_unphased(1);
                        }
                        else{
                            species_gt[j * 2] = bcf_gt_unphased(1);
                            species_gt[j * 2 + 1] = bcf_gt_unphased(1);
                        }
                    }
                }
                // Insert new genotypes into output record
                bcf_update_genotypes(out_hdr, out_rec, species_gt.data(), species_gt.size());
                int ret = bcf_write1(outf, out_hdr, out_rec);
                bcf_clear(out_rec);
                kept++;
            }
        }
        ++nsnp;
        if (nsnp % 1000 == 0){
            fprintf(stderr, "Processed %ld SNPs\r", nsnp);
        }
    }
    
    hts_close(bcf_reader);
    bcf_destroy(bcf_record);
    bcf_destroy(out_rec);
    bcf_hdr_destroy(bcf_header);
    bcf_hdr_destroy(out_hdr);
    hts_close(outf);
    
    fprintf(stderr, "Processed %ld SNPs\n", nsnp);
    fprintf(stderr, "Kept %ld of %ld SNPs (%.2f%%)\n", kept, nsnp, 100.0*(double)kept/(double)nsnp);

    return 0;
}
