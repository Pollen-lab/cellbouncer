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
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdlib>
#include <htslib/kseq.h>
#include <htswrapper/bc.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "reads_demux.h"

using namespace std;

KSEQ_INIT(gzFile, gzread);

// Constructor
reads_demuxer::reads_demuxer(bc_whitelist& whitelist,
    robin_hood::unordered_map<unsigned long, short>& bc2species,
    map<short, string>& idx2species,
    string& outdir){
    
    this->whitelist = whitelist;
    this->bc2species = bc2species;
    this->idx2species = idx2species;
    this->outdir = outdir;

    initialized = false;
}

// Destructor
reads_demuxer::~reads_demuxer(){
    if (initialized){
        for (int i = 0; i < outfiles.size(); ++i){
            gzclose(outfiles[i]);
        }
    }
}

// If already initialized for a given set of files, close.
void reads_demuxer::close(){
    if (this->initialized){
        for (int i = 0; i < outfiles.size(); ++i){
            gzclose(outfiles[i]);
        }
        r1 = "";
        r2 = "";
        r3 = "";
        is_atac = false;
         outfiles.clear(); 
        initialized = false;
    }
}

// Set up for RNA-seq reads or feature barcoding data
void reads_demuxer::init_rna_or_custom(string file_prefix, string& r1filename, string& r2filename){
    
    // If already initialized for something else, close existing files
    close();

    // Trim paths from filenames for output filenames
    string r1filetrim = filename_nopath(r1filename);
    string r2filetrim = filename_nopath(r2filename);

    // Append prefix (i.e. "GEX") identifier to output files, if necessary
    string prefix = file_prefix + "_";
    if (r1filetrim.length() < prefix.length() || r1filetrim.substr(0, prefix.length()) != prefix){
        r1filetrim = prefix + r1filetrim;
    }
    if (r2filetrim.length() < prefix.length() || r2filetrim.substr(0, prefix.length()) != prefix){
        r2filetrim = prefix + r2filetrim;
    }

    // To keep the 10X pipeline happy, we'll create output files with the same name as
    // input files, separated in directories according to species of origin.
    if (outdir[outdir.length()-1] != '/'){
        outdir += "/";
    }
    
    outfiles.reserve(2 * idx2species.size());

    for (map<short, string>::iterator spec = idx2species.begin(); spec != idx2species.end(); ++spec){
        string dirn = outdir + spec->second;
        if (!mkdir(dirn.c_str(), 0775)){
            // Assume directory already exists
        }
        string fn = dirn + "/" + r1filetrim;
        outfiles[spec->first * 2] = gzopen(fn.c_str(), "w");
        if (!outfiles[spec->first * 2]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
        fn = dirn + "/" + r2filetrim;
        outfiles[spec->first * 2 + 1] = gzopen(fn.c_str(), "w");
        if (!outfiles[spec->first * 2 + 1]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
    }

    r1 = r1filename;
    r2 = r2filename;
    r3 = "";
    initialized = true;
}

void reads_demuxer::init_rna(string& r1filename, string& r2filename){
    init_rna_or_custom("GEX", r1filename, r2filename);
}

void reads_demuxer::init_custom(string prefix, string& r1filename, string& r2filename){
    init_rna_or_custom(prefix, r1filename, r2filename);
}

bool reads_demuxer::scan_rna(){
    if (!initialized || is_atac){
        return false;
    }
    
    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    gzFile f_fp;
    gzFile r_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    f_fp = gzopen(r1.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1.c_str());
        exit(1);
    }    
    r_fp = gzopen(r2.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2.c_str());
        exit(1);
    }
   
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    
    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R2 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        unsigned long bc_key;
        if (whitelist.lookup1_bf(seq_f->seq.s, bc_key, false)){
            int species = bc2species[bc_key];
            write_fastq(seq_f->name.s, seq_f->name.l, seq_f->seq.s, seq_f->seq.l, seq_f->qual.s,
                species*2);
            write_fastq(seq_f->name.s, seq_r->name.l, seq_r->seq.s, seq_r->seq.l, seq_r->qual.s,
                species*2 + 1);
        } 
    }

    kseq_destroy(seq_f);
    kseq_destroy(seq_r);
    return true;  
}

bool reads_demuxer::scan_custom(){
    return scan_rna();
}

/**
 * Set up for ATAC-seq files
 */
void reads_demuxer::init_atac(string& r1filename, string& r2filename, string& r3filename){
    close();
    
    // Trim directory paths off file names for output file names
    string r1filetrim = filename_nopath(r1filename);
    string r2filetrim = filename_nopath(r2filename);
    string r3filetrim = filename_nopath(r3filename);
    
    // If no "ATAC" identifier present in filenames, add it
    if (r1filetrim.length() < 5 || r1filetrim.substr(0, 5) != "ATAC_"){
        r1filetrim = "ATAC_" + r1filetrim;
    }
    if (r2filetrim.length() < 5 || r2filetrim.substr(0, 5) != "ATAC_"){
        r2filetrim = "ATAC_" + r2filetrim;
    }
    if (r3filetrim.length() < 5 || r3filetrim.substr(0, 5) != "ATAC_"){
        r3filetrim = "ATAC_" + r3filetrim;
    }

    // To keep the 10X pipeline happy, we'll create output files with the same name as
    // input files, separated in directories according to species of origin.
    if (outdir[outdir.length()-1] != '/'){
        outdir += "/";
    }
    
    outfiles.reserve(idx2species.size() * 3);

    for (map<short, string>::iterator spec = idx2species.begin(); spec != idx2species.end(); ++spec){
        string dirn = outdir + spec->second;
        if (!mkdir(dirn.c_str(), 0775)){
            // Assume directory already exists
        }
        string fn = dirn + "/" + r1filetrim;
        outfiles[spec->first * 3] = gzopen(fn.c_str(), "w");
        if (!outfiles[spec->first * 3]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
        fn = dirn + "/" + r2filetrim;
        outfiles[spec->first * 3 + 1] = gzopen(fn.c_str(), "w");
        if (!outfiles[spec->first * 3 + 1]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
        fn = dirn + "/" + r3filetrim;
        outfiles[spec->first * 3 + 2] = gzopen(fn.c_str(), "w");
        if (!outfiles[spec->first * 3 + 2]){
            fprintf(stderr, "ERROR opening %s for writing.\n", fn.c_str());
            exit(1);
        }
    }
    r1 = r1filename;
    r2 = r2filename;
    r3 = r3filename;
    initialized = true;
}

bool reads_demuxer::scan_atac(){
    if (!initialized || !is_atac){
        return false;
    }

    // Now iterate through read files, find/match barcodes, and assign to the correct files.
    // Prep input file(s).
    int f_progress;
    int r_progress;
    int i_progress;
    gzFile f_fp;
    gzFile r_fp;
    gzFile i_fp;
    kseq_t* seq_f;
    kseq_t* seq_r;
    kseq_t* seq_i;
    f_fp = gzopen(r1.c_str(), "r");
    if (!f_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r1.c_str());
        exit(1);
    }    
    r_fp = gzopen(r3.c_str(), "r");
    if (!r_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r3.c_str());
        exit(1);
    }
    i_fp = gzopen(r2.c_str(), "r");
    if (!i_fp){
        fprintf(stderr, "ERROR opening %s for reading\n", r2.c_str());
        exit(1);
    }
    seq_f = kseq_init(f_fp);
    seq_r = kseq_init(r_fp);
    seq_i = kseq_init(i_fp);

    while ((f_progress = kseq_read(seq_f)) >= 0){
        r_progress = kseq_read(seq_r);
        if (r_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R3 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        i_progress = kseq_read(seq_i);
        if (i_progress < 0){
            fprintf(stderr, "ERROR: read order not matching R2 at seq %s in R1 file\n", seq_f->name.s);
            exit(1);
        }
        
        unsigned long bc_key;
        // ATAC barcode is in seq_i
        // Orientation: reverse complement, end of read
        if (whitelist.lookup2_er(seq_i->seq.s, bc_key, false)){
            int species = bc2species[bc_key];
            write_fastq(seq_f->name.s, seq_f->name.l, seq_f->seq.s, 
                seq_f->seq.l, seq_f->qual.s, species*3);
            write_fastq(seq_i->name.s, seq_i->name.l, seq_i->seq.s, 
                seq_i->seq.l, seq_i->qual.s, species*3+1);
            write_fastq(seq_r->name.s, seq_r->name.l, seq_r->seq.s, 
                seq_r->seq.l, seq_r->qual.s, species*3+2);
        }
    } 

    kseq_destroy(seq_f);
    kseq_destroy(seq_r);
    kseq_destroy(seq_i);
    return true;
}

/**
 * Write a FASTQ record to disk.
 */
void reads_demuxer::write_fastq(const char* id, 
    int idlen, 
    const char* seq, 
    int seqlen, 
    const char* qual, 
    int out_idx){
    
    gzwrite(outfiles[out_idx], "@", 1);
    gzwrite(outfiles[out_idx], id, idlen);
    gzwrite(outfiles[out_idx], "\n", 1);
    gzwrite(outfiles[out_idx], seq, seqlen);
    gzwrite(outfiles[out_idx], "\n", 1);
    gzwrite(outfiles[out_idx], "+\n", 2);
    gzwrite(outfiles[out_idx], qual, seqlen);
    gzwrite(outfiles[out_idx], "\n", 1);
}
