#include <zlib.h>
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
#include <unordered_map>
#include <set>
#include <deque>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <htslib/kseq.h>
#include <zlib.h>
#include <limits.h>
#include <thread>
#include <sys/stat.h>
#include <condition_variable>
#include <mutex>
#include <htswrapper/bc.h>
#include <htswrapper/umi.h>
#include <htswrapper/gzreader.h>
#include <htswrapper/khashtable.h>
#include <random>
#include "common.h"
#include "species_kmers.h"

KSEQ_INIT(gzFile, gzread);

using std::cout;
using std::endl;
using namespace std;

// ===== species_kmers.cpp
// This file contains functions used to scan reads for species-specific k-mers
// Used by demux_species.

rp_info::~rp_info(){
    //free(seq_f);
    //free(seq_r);  
}

rp_info::rp_info(const char* kseq_f, int seq_f_len, const char* kseq_r, int seq_r_len){
    this->seq_f = (char*)malloc((seq_f_len+1)*sizeof(char));
    this->seq_r = (char*)malloc((seq_r_len+1)*sizeof(char));
    strncpy(this->seq_f, kseq_f, seq_f_len);
    strncpy(this->seq_r, kseq_r, seq_r_len);
    this->seq_f[seq_f_len] = '\0';
    this->seq_r[seq_r_len] = '\0';
    this->len_f = seq_f_len;
    this->len_r = seq_r_len;
}

rp_info::rp_info(const rp_info& other){
    this->seq_f = other.seq_f;
    this->seq_r = other.seq_r;
    this->len_f = other.len_f;
    this->len_r = other.len_r;
}

/*
kmer_node_ptr::kmer_node_ptr(){
    f_A = NULL;
    f_C = NULL;
    f_G = NULL;
    f_T = NULL;
    f_A_flip = false;
    f_C_flip = false;
    f_G_flip = false;
    f_T_flip = false;
    r_A = NULL;
    r_C = NULL;
    r_G = NULL;
    r_T = NULL;
    r_A_flip = false;
    r_C_flip = false;
    r_G_flip = false;
    r_T_flip = false;
}
*/
species_kmer_counter::species_kmer_counter(int nt, 
    int k, 
    int ns,
    bc_whitelist* wl,
    robin_hood::unordered_map<unsigned long, map<short, int> >* bsc,
    int umi_start,
    int umi_len): tab(k){
    
    this->terminate_threads = false;
    this->num_threads = nt;
    this->k = k;
    this->wl = wl;
    this->num_species = ns;
    this->initialized = false;
    this->on = false;
    this->bc_species_counts = bsc;
    this->umi_start = umi_start;
    this->umi_len = umi_len;
    
    this->n_samp = -1;

    this->use_umis = true;
    if (umi_start == -1 || umi_len == -1){
        use_umis = false;
    }
}

species_kmer_counter::~species_kmer_counter(){
    if (this->on){
        this->close_pool();
        this->on = false;
    }
    if (this->initialized){
        //kmsuftree_destruct(kt, 0);
        tab.clear();
    }
    // Free all UMI counters
    for (robin_hood::unordered_map<unsigned long, umi_set_exact* >::iterator x = bc_species_umis.begin();
        x != bc_species_umis.end(); ++x){
        delete x->second;
        x->second = NULL;
    }
    bc_species_umis.clear();
}

void species_kmer_counter::set_n_samp(int ns){
    this->n_samp = ns;
}

void species_kmer_counter::init(short species_idx, string& kmerfile){
    this->terminate_threads = false;

    //this_species = species_idx;

    if (this->initialized){
        tab.clear();
        //kmsuftree_destruct(kt, 0);
    }
    //tab = khashtable<short>(k); 
    parse_kmer_counts_serial(kmerfile, species_idx);
    
    this->initialized = true;
}

void species_kmer_counter::add(short species_idx, string& kmerfile){
    if (!this->initialized){
        //tab = khashtable<short>(k);
        //kt = init_kmer_tree(k);
        this->initialized = true;
        this->terminate_threads = false;
    }
    parse_kmer_counts_serial(kmerfile, species_idx);
}

/**
 * Only using canonical k-mers = first in sort order relative to reverse complement.
 * Check if a given k-mer or its reverse complement comes first in alphabetical sort order.
 *
bool species_kmer_counter::is_rc_first(const char* kmer, int k){
    for (int i = 0; i < k; ++i){
        char char_f = kmer[k];
        switch(kmer[k]){
            case 'a':
            case 'A':
                char_f = 'A';
                break;
            case 'c':
            case 'C':
                char_f = 'C';
                break;
            case 'g':
            case 'G':
                char_f = 'G';
                break;
            case 't':
            case 'T':
                char_f = 'T';
                break;
            default:
                char_f = 'N';
                break;
        }
        char char_rc = 'N';
        switch(kmer[k-1-i]){
            case 'A':
            case 'a':
                char_rc = 'T';
                break;
            case 'c':
            case 'C':
                char_rc = 'G';
                break;
            case 'g':
            case 'G':
                char_rc = 'C';
                break;
            case 't':
            case 'T':
                char_rc = 'A';
                break;
            default:
                char_rc = 'N';
                break;
        }
        if (char_f == 'N' || char_rc == 'N'){
            // No matter; neither will be in the lookup table.
            return false;
        }
        if (char_rc < char_f){
            return true;
        }
        else if (char_f < char_rc){
            return false;
        }   
    }
    return false;
}
*/

void species_kmer_counter::disable_umis(){
    this->use_umis = false;
}

void species_kmer_counter::enable_umis(){
    this->use_umis = true;
}

/**
 * Adds sequences from reads to data that will be retrieved by worker threads.
 */
void species_kmer_counter::process_gex_files(string& r1filename, 
    string& r2filename){
    
    if (num_threads > 1){
        launch_gex_threads();
    }
    else{
        vector<int> v;
        for (int i = 0; i < num_species; ++i){
            v.push_back(0);
        }
        species_counts.push_back(v);
        khashkeys.emplace_back(k);
    }

    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    gzFile f_fp;
    gzFile r_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    f_fp = gzopen(r1filename.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1filename.c_str());
        exit(1);
    }    
    r_fp = gzopen(r2filename.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2filename.c_str());
        exit(1);
    }
   
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: %s still contains reads, but %s reached end of file\n", 
                r1filename.c_str(), r2filename.c_str());
            fprintf(stderr, "%s is likely truncated or corrupted.\n", r2filename.c_str());
            exit(1);
        }
        if (num_threads > 1){
            add_rp_job(seq_f->seq.s, seq_f->seq.l, seq_r->seq.s, seq_r->seq.l);        
        }
        else{
            // Just count normally, without wasting overhead counting sequences
            scan_gex_data(seq_f->seq.s, seq_f->seq.l, seq_r->seq.s, seq_r->seq.l, 0);
        }
    }

    // Check that R2 has hit EOF as well.
    if ((r_progress = kseq_read(seq_r)) >= 0){
        fprintf(stderr, "ERROR: %s still contains reads, but %s reached end of file\n", 
            r2filename.c_str(), r1filename.c_str());
        fprintf(stderr, "%s is likely truncated or corrupted.\n", r1filename.c_str());
        exit(1);
    }
    
    if (num_threads > 1){
       close_pool();
    }

    for (robin_hood::unordered_map<unsigned long, umi_set_exact* >::iterator x = bc_species_umis.begin();
        x != bc_species_umis.end(); ++x){
        delete x->second;
        x->second = NULL;
    }
    bc_species_umis.clear();
}

// Stop all running threads and then destroy them
void species_kmer_counter::close_pool(){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->terminate_threads = true;   
    }
    this->has_jobs.notify_all();
    for (int i = 0; i < this->threads.size(); ++i){
        this->threads[i].join();
    }
    this->threads.clear();
    this->on = false;
}

/**
 * Add to the queue of jobs to read paired reads
 */
void species_kmer_counter::add_rp_job(const char* seq_f, int seq_f_len, const char* seq_r, int seq_r_len){

    {
        // Add job to processing queue
        unique_lock<mutex> lock(this->queue_mutex);
        this->rp_jobs.emplace_back(seq_f, seq_f_len, seq_r, seq_r_len);
    }
    this->has_jobs.notify_one();
}

/**
 * Add to the queue of jobs to read read triplets
 */
void species_kmer_counter::add_rt_job(rt_info& info){
    {
        unique_lock<mutex> lock(this->queue_mutex);
        this->rt_jobs.push_back(info);
    }
    this->has_jobs.notify_one();
}

// Count k-mers for one species in a specific read.
void species_kmer_counter::scan_seq_kmers(const char* seq, int len, int* result_counts, khashkey& key){
    key.reset();
    int pos = 0;
    while(key.scan_kmers(seq, len, pos)){
        short spec;
        if (tab.lookup(key, spec)){
            result_counts[spec]++;
            return;
        }
        else{
            //key.print(false);
        }
    }
    /*
    kmer_node_ptr* cur = NULL;
    bool cur_rc = false;

    static int tot = 0;
    static int pointed = 0;
    static int lookedup = 0;
    
    int N_pos = -1;
    
    size_t k_idx;
    size_t k_idx_rev;
    bool rebuild = true;

    for (int i = 0; i < len - kt->k; ++i){
            
        // If contains N, then can't have valid k-mer.
        if (seq[i + kt->k - 1] == 'N'){
            N_pos = i + kt->k - 1;
        }
        if (N_pos >= i){
            // N will occur in this k-mer.
            cur = NULL;
        }
        else{
            if (cur == NULL){
                // Attempt to look up.
                void* dat = NULL;
                if (rebuild){
                    if ( kmer2inx( seq + i, kt->k_ar_size, &k_idx, 0, kt->k ) ){
                        kmer2inx( seq + i, kt->k_ar_size, &k_idx_rev, 1, kt->k );
                        rebuild = false;
                        if (is_rc_first( seq + i, kt->k )){
                            dat = kmer_tree_lookup_aux(seq + i, kt, 1, &k_idx_rev);
                            if (dat != NULL){
                                cur_rc = true;
                            }
                        }
                        else{
                            dat = kmer_tree_lookup_aux(seq + i, kt, 0, &k_idx);  
                            if (dat != NULL){
                                cur_rc = false;
                            }
                        }
                    }
                    else{
                        // Still need to rebuild.
                    }
                }
                else{
                    if ( kmer2inx_continue( seq + i, kt->k_ar_size, &k_idx, 0, kt->k ) ){
                        kmer2inx_continue( seq + i, kt->k_ar_size, &k_idx_rev, 1, kt->k );
                        rebuild = false;
                        if (is_rc_first( seq + i, kt->k )){
                            dat = kmer_tree_lookup_aux(seq + i, kt, 1, &k_idx_rev);
                            if (dat != NULL){
                                cur_rc = true;
                            }
                        }
                        else{
                            dat = kmer_tree_lookup_aux(seq + i, kt, 0, &k_idx);
                            if (dat != NULL){
                                cur_rc = false;
                            }
                        }
                    }
                    else{
                        // N encountered; rebuild
                        rebuild = true;
                    }
                }        
                if (dat != NULL){
                    //++species_count;
                    cur = (kmer_node_ptr*)dat;
                    result_counts[cur->species]++;
                    lookedup++;
                    return;
               } 
            }
            else{
                // We already have a k-mer to start from. See if we can follow
                // the path.
                char newbase = seq[i + kt->k - 1];
                if (cur_rc){
                    if (newbase == 'A'){
                        if (cur->r_A != NULL){
                            cur_rc = !cur->r_A_flip;
                            cur = cur->r_A;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                    else if (newbase == 'C'){
                        if (cur->r_C != NULL){
                            cur_rc = !cur->r_C_flip;
                            cur = cur->r_C;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                    else if (newbase == 'G'){
                        if (cur->r_G != NULL){
                            cur_rc = !cur->r_G_flip;
                            cur = cur->r_G;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{

                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                    else if (newbase == 'T'){
                        if (cur->r_T != NULL){
                            cur_rc = !cur->r_T_flip;
                            cur = cur->r_T;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                }
                else{
                    if (newbase == 'A'){
                        if (cur->f_A != NULL){
                            cur = cur->f_A;        
                            cur_rc = cur->f_A_flip;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                    else if (newbase == 'C'){
                        if (cur->f_C != NULL){
                            cur = cur->f_C;
                            cur_rc = cur->f_C_flip;
                            result_counts[cur->species]++;
                            //species_count++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                    else if (newbase == 'G'){
                        if (cur->f_G != NULL){
                            cur = cur->f_G;
                            cur_rc = cur->f_G_flip;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                    else if (newbase == 'T'){
                        if (cur->f_T != NULL){
                            cur = cur->f_T;
                            cur_rc = cur->f_T_flip;
                            //species_count++;
                            result_counts[cur->species]++;
                            pointed++;
                        }
                        else{
                            cur = NULL;
                            cur_rc = false;
                        }
                    }
                }
            }
        }
    }
    */
}

void species_kmer_counter::scan_gex_data(const char* seq_f, 
    int seq_f_len, const char* seq_r, int seq_r_len, int thread_idx){
    
    unsigned long bc_key = 0;
    bool exact;
    
    if (wl->lookup(seq_f, bc_key, exact, seq_f_len)){
        
        bool dup_read = false;
        if (use_umis && umi_start >= 0 && umi_len > 0){
            umi u(seq_f + umi_start, umi_len);
            
            umi_set_exact* us; 
            {
               unique_lock<mutex> lock(umi_mutex);
                 
               if (bc_species_umis.count(bc_key) == 0){
                   us = new umi_set_exact();
                   bc_species_umis.emplace(bc_key, us); 
                }
                else{
                    us = bc_species_umis[bc_key];
                }
            }
            {
                unique_lock<mutex> lock(us->umi_mutex);
                dup_read = us->add(u);
            }
        }
        
        if (dup_read){
            return;
        }

        // In 10x scRNA-seq data, only the reverse read contains information
        // forward read is barcode
        for (int j = 0; j < num_species; ++j){
            species_counts[thread_idx][j] = 0;
        }
        scan_seq_kmers(seq_r, seq_r_len, species_counts[thread_idx].data(), khashkeys[thread_idx]);
        
        for (int j = 0; j < num_species; ++j){
            if (species_counts[thread_idx][j] > 0){
                int this_species = j;
                int nk = species_counts[thread_idx][j];
                unique_lock<mutex> lock(this->counts_mutex);
                if (bc_species_counts->count(bc_key) == 0){
                    map<short, int> m;
                    bc_species_counts->emplace(bc_key, m);
                }
                if ((*bc_species_counts)[bc_key].count(this_species) == 0){
                    (*bc_species_counts)[bc_key].insert(make_pair(this_species, 0));
                }
                //(*bc_species_counts)[bc_key][this_species]++;
                (*bc_species_counts)[bc_key][this_species] += nk;
            }
        }
    }
}

char complement(char base){
    if (base == 'A'){
        return 'T';
    }
    else if (base == 'C'){
        return 'G';
    }
    else if (base == 'G'){
        return 'C';
    }
    else if (base == 'T'){
        return 'A';
    }
    return 'N';
}

/**
 * Assumes we're adding both forward and RC versions of every k-mer.
 * User must reverse complement the sequence and then call this function;
 * assumes forward orientation.
 *
void species_kmer_counter::add_node_to_tree_simple(char* kmer, char* kmer_buf, kmer_tree_p kt, short species_idx){
    char bases[4] = {'A', 'C', 'G', 'T'};

    kmer_node_ptr* n = (kmer_node_ptr*) malloc(sizeof(kmer_node_ptr));
    // Initialize.
    n->f_A = NULL;
    n->f_C = NULL;
    n->f_G = NULL;
    n->f_T = NULL;
    n->r_A = NULL;
    n->r_C = NULL;
    n->r_G = NULL;
    n->r_T = NULL;
    n->f_A_flip = false;
    n->f_C_flip = false;
    n->f_G_flip = false;
    n->f_T_flip = false;
    n->r_A_flip = false;
    n->r_C_flip = false;
    n->r_G_flip = false;
    n->r_T_flip = false;
    n->species = species_idx;

    // Check whether any 1-base-away k-mer (moving backward or forward)
    //  has already been created. Link them up if so.

    // Rules: moving forward in sequence:
    //  b1 = last base being added
    //  b2 = complement of first base
    //  is current k-mer rc? 
    //      yes -> link from cur k-mer b1, reverse
    //          flip = !next_rc
    //      no -> link from cur k-mer b1, forward
    //          flip = next_rc
    //  is next k-mer rc?
    //      yes -> link from next k-mer b2, forward
    //          flip = !cur_rc
    //      no -> link from next k-mer b2, reverse
    //          flip = cur_rc
    
    // Check moving forward
    strncpy(&kmer_buf[0], kmer + 1, kt->k-1);
    kmer_buf[kt->k] = '\0';
    
    for (int i = 0; i < 4; ++i){
        kmer_node_ptr* ptr = NULL; 
    
        kmer_buf[kt->k-1] = bases[i];
        
        void* res = kmer_tree_lookup(&kmer_buf[0], kt, 0);
        if (res != NULL){
            ptr = (kmer_node_ptr*)res;
        }
        
        if (ptr != NULL){
            if (bases[i] == 'A'){
                n->f_A = ptr;
            }
            else if (bases[i] == 'C'){
                n->f_C = ptr;
            }
            else if (bases[i] == 'G'){
                n->f_G = ptr;
            }
            else if (bases[i] == 'T'){
                n->f_T = ptr;
            }
            if (kmer[0] == 'A'){
                ptr->r_T = n;
            }
            else if (kmer[0] == 'C'){
                ptr->r_G = n;
            }
            else if (kmer[0] == 'G'){
                ptr->r_C = n;
            }
            else if (kmer[0] == 'T'){
                ptr->r_A = n;
            }
        }
    }

    // Check moving backward
    strncpy(&kmer_buf[1], kmer, kt->k-1);
    kmer_buf[kt->k] = '\0';
    for (int i = 0; i < 4; ++i){
        kmer_node_ptr* ptr = NULL; 
        kmer_buf[0] = bases[i];
        
        void* res = kmer_tree_lookup(&kmer_buf[0], kt, 0);
        if (res != NULL){
            ptr = (kmer_node_ptr*)res;
        }
        
        if (ptr != NULL){
            // Same rules as above, but now buf and line are reversed.
            if (kmer[kt->k-1] == 'A'){
                ptr->f_A = n;
            }
            else if (kmer[kt->k-1] == 'C'){
                ptr->f_C = n;
            }
            else if (kmer[kt->k-1] == 'G'){
                ptr->f_G = n;
            }
            else if (kmer[kt->k-1] == 'T'){
                ptr->f_T = n;
            }
            if (bases[i] == 'A'){
                n->r_T = ptr;
            }
            else if (bases[i] == 'C'){
                n->r_G = ptr;
            }
            else if (bases[i] == 'G'){
                n->r_C = ptr;
            }
            else if (bases[i] == 'T'){
                n->r_A = ptr;
            }
        }
    }

    // Store the node in the suffix array.
    kmer_tree_add(kmer, kt, (void*)n, false);
}
*/
/** 
 * Only adds nodes in the orientation (F,RC) that comes alphabetically first
 *
void species_kmer_counter::add_node_to_tree_canonical(char* kmer, char* buf, kmer_tree_p kt, short species_idx){
    char bases[4] = {'A', 'C', 'G', 'T'};

    kmer_node_ptr* n = (kmer_node_ptr*) malloc(sizeof(kmer_node_ptr));
    // Initialize.
    n->f_A = NULL;
    n->f_C = NULL;
    n->f_G = NULL;
    n->f_T = NULL;
    n->r_A = NULL;
    n->r_C = NULL;
    n->r_G = NULL;
    n->r_T = NULL;
    n->f_A_flip = false;
    n->f_C_flip = false;
    n->f_G_flip = false;
    n->f_T_flip = false;
    n->r_A_flip = false;
    n->r_C_flip = false;
    n->r_G_flip = false;
    n->r_T_flip = false;
    n->species = species_idx;

    bool rc = is_rc_first(kmer, kt->k);
    
    // Check whether any 1-base-away k-mer (moving backward or forward)
    //  has already been created. Link them up if so.

    // Rules: moving forward in sequence:
    //  b1 = last base being added
    //  b2 = complement of first base
    //  is current k-mer rc? 
    //      yes -> link from cur k-mer b1, reverse
    //          flip = !next_rc
    //      no -> link from cur k-mer b1, forward
    //          flip = next_rc
    //  is next k-mer rc?
    //      yes -> link from next k-mer b2, forward
    //          flip = !cur_rc
    //      no -> link from next k-mer b2, reverse
    //          flip = cur_rc

    // Check moving forward
    strncpy(&buf[0], kmer + 1, kt->k-1);
    buf[kt->k] = '\0';
    
    for (int i = 0; i < 4; ++i){
        kmer_node_ptr* ptr = NULL; 
        bool ptr_rc = false;
    
        buf[kt->k-1] = bases[i];
        
        if (is_rc_first(&buf[0], kt->k)){
            void* res = kmer_tree_lookup(&buf[0], kt, 1);
            if (res != NULL){
                ptr = (kmer_node_ptr*)res;
                ptr_rc = true;
            }
        }
        else{
            void* res = kmer_tree_lookup(&buf[0], kt, 0);
            if (res != NULL){
                ptr = (kmer_node_ptr*)res;
                ptr_rc = false;
            }
        }
        if (ptr != NULL){
            if (rc){
                if (bases[i] == 'A'){
                    n->r_A = ptr;
                    n->r_A_flip = !ptr_rc;
                }   
                else if (bases[i] == 'C'){
                    n->r_C = ptr;
                    n->r_C_flip = !ptr_rc;
                }
                else if (bases[i] == 'G'){
                    n->r_G = ptr;
                    n->r_G_flip = !ptr_rc;
                }
                else if (bases[i] == 'T'){
                    n->r_T = ptr;
                    n->r_T_flip = !ptr_rc;
                }
            }
            else{
                if (bases[i] == 'A'){
                    n->f_A = ptr;
                    n->f_A_flip = ptr_rc;
                }
                else if (bases[i] == 'C'){
                    n->f_C = ptr;
                    n->f_C_flip = ptr_rc;
                }
                else if (bases[i] == 'G'){
                    n->f_G = ptr;
                    n->f_G_flip = ptr_rc;
                }
                else if (bases[i] == 'T'){
                    n->f_T = ptr;
                    n->f_T_flip = ptr_rc;
                }
            }
            if (ptr_rc){
                if (kmer[0] == 'A'){
                    ptr->f_T = n;
                    ptr->f_T_flip = !rc;
                }
                else if (kmer[0] == 'C'){
                    ptr->f_G = n;
                    ptr->f_G_flip = !rc;
                }
                else if (kmer[0] == 'G'){
                    ptr->f_C = n;
                    ptr->f_C_flip = !rc;
                }
                else if (kmer[0] == 'T'){
                    ptr->f_A = n;
                    ptr->f_A_flip = !rc;
                }
            }
            else{
                if (kmer[0] == 'A'){
                    ptr->r_T = n;
                    ptr->r_T_flip = rc;
                }
                else if (kmer[0] == 'C'){
                    ptr->r_G = n;
                    ptr->r_G_flip = rc;
                }
                else if (kmer[0] == 'G'){
                    ptr->r_C = n;
                    ptr->r_C_flip = rc;
                }
                else if (kmer[0] == 'T'){
                    ptr->r_A = n;
                    ptr->r_A_flip = rc;
                }
            }
        }
    }

    // Check moving backward
    strncpy(&buf[1], kmer, kt->k-1);
    buf[kt->k] = '\0';
    for (int i = 0; i < 4; ++i){
        kmer_node_ptr* ptr = NULL; 
        bool ptr_rc = false;
        buf[0] = bases[i];
        
        if (is_rc_first(&buf[0], kt->k)){
            void* res = kmer_tree_lookup(&buf[0], kt, 1);
            if (res != NULL){
                ptr = (kmer_node_ptr*)res;
                ptr_rc = true;
            }
        }
        else{
            void* res = kmer_tree_lookup(&buf[0], kt, 0);
            if (res != NULL){
                ptr = (kmer_node_ptr*)res;
                ptr_rc = false;
            }
        }
        if (ptr != NULL){
            // Same rules as above, but now buf and line are reversed.
            if (ptr_rc){
                if (kmer[k-1] == 'A'){
                    ptr->r_A = n;
                    ptr->r_A_flip = !rc;
                }   
                else if (kmer[k-1] == 'C'){
                    ptr->r_C = n;
                    ptr->r_C_flip = !rc;
                }
                else if (kmer[k-1] == 'G'){
                    ptr->r_G = n;
                    ptr->r_G_flip = !rc;
                }
                else if (kmer[k-1] == 'T'){
                    ptr->r_T = n;
                    ptr->r_T_flip = !rc;
                }
            }
            else{
                if (kmer[k-1] == 'A'){
                    ptr->f_A = n;
                    ptr->f_A_flip = rc;
                }
                else if (kmer[k-1] == 'C'){
                    ptr->f_C = n;
                    ptr->f_C_flip = rc;
                }
                else if (kmer[k-1] == 'G'){
                    ptr->f_G = n;
                    ptr->f_G_flip = rc;
                }
                else if (kmer[k-1] == 'T'){
                    ptr->f_T = n;
                    ptr->f_T_flip = rc;
                }
            }
            if (rc){
                if (bases[i] == 'A'){
                    n->f_T = ptr;
                    n->f_T_flip = !ptr_rc;
                }
                else if (bases[i] == 'C'){
                    n->f_G = ptr;
                    n->f_G_flip = !ptr_rc;
                }
                else if (bases[i] == 'G'){
                    n->f_C = ptr;
                    n->f_C_flip = !ptr_rc;
                }
                else if (bases[i] == 'T'){
                    n->f_A = ptr;
                    n->f_A_flip = !ptr_rc;
                }
            }
            else{
                if (bases[i] == 'A'){
                    n->r_T = ptr;
                    n->r_T_flip = ptr_rc;
                }
                else if (bases[i] == 'C'){
                    n->r_G = ptr;
                    n->r_G_flip = ptr_rc;
                }
                else if (bases[i] == 'G'){
                    n->r_C = ptr;
                    n->r_C_flip = ptr_rc;
                }
                else if (bases[i] == 'T'){
                    n->r_A = ptr;
                    n->r_A_flip = ptr_rc;
                }
            }
        }
    }

    // Store the node in the suffix array.
    kmer_tree_add(kmer, kt, (void*)n, rc);
}
*/

void species_kmer_counter::parse_kmer_counts_serial(string& countsfilename,
    short species_idx){
    
    gzreader reader(countsfilename);
    int k = -1;
    
    while (reader.next()){
        if (k == -1){
            k = strlen(reader.line);
        }
        tab.add(reader.line, species_idx);
        //add_node_to_tree_canonical(reader.line, buf2, kt, species_idx);
    }
}

/**
 * Worker function for a GEX read scanning job
 */
void species_kmer_counter::gex_thread(int thread_idx){
    
     while(true){
        char* seq_f = NULL;
        int seq_f_len = 0;
        char* seq_r = NULL;
        int seq_r_len = 0;
        {
            unique_lock<mutex> lock(this->queue_mutex);
            this->has_jobs.wait(lock, [this]{ return rp_jobs.size() > 0 ||
                terminate_threads;});
            if (this->rp_jobs.size() == 0 && this->terminate_threads){
                return;
            }
            seq_f = this->rp_jobs[0].seq_f;
            seq_f_len = this->rp_jobs[0].len_f;
            seq_r = this->rp_jobs[0].seq_r;
            seq_r_len = this->rp_jobs[0].len_r;
            
            // Can't free here; need pointers to stick around
            this->rp_jobs.pop_front();
        }
        if (seq_f != NULL && seq_r != NULL){
            scan_gex_data(seq_f, seq_f_len, seq_r, seq_r_len, thread_idx);
            free(seq_f);
            free(seq_r);
        }
     }
 }

void species_kmer_counter::launch_gex_threads(){
    
    if (!this->initialized){
        fprintf(stderr, "ERROR: not initialized\n");
        exit(1);
    }

    this->terminate_threads = false;

    for (int i = 0; i < this->num_threads; ++i){
        
        vector<int> v;
        this->species_counts.push_back(v);
        for (int j = 0; j < num_species; ++j){
            this->species_counts[this->species_counts.size()-1].push_back(0);
        }
        this->khashkeys.emplace_back(k);

        this->threads.push_back(thread(&species_kmer_counter::gex_thread,
            this, i));
        
    
    }

    this->on = true;
}


