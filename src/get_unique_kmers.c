#include <getopt.h>
#include <argp.h>
#include <zlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "libfastk.h"

KSEQ_INIT(gzFile, gzread);

void help(int code){
    fprintf(stderr, "USAGE: get_unique_kmers\n");
    fprintf(stderr, "After running FASTK to count k-mers on different species' transcriptomes, \
run this prgram to read FASTK's k-mer table files and output lists of k-mers unique to each \
species. These can then be used to demultiplex reads from multiple species without mapping \
to a reference genome.\n");
    fprintf(stderr, "OPTIONS:\n");
    fprintf(stderr, "--kmers -k Two or more FASTK tables from reference genomes (specify \
-k multiple times)\n");
    fprintf(stderr, "--output_prefix -o Where to write species-specific kmers. Files will \
be in the format {outprefix}.{index}.kmers.\n"); 
    exit(code);
}

int kcomp(const Kmer_Stream* k1, const Kmer_Stream* k2){
    if (k1->cpre < k2->cpre){
        return -1;
    }
    else if (k1->cpre > k2->cpre){
        return 1;
    }
    else{
        for (int idx = 0; idx < k1->kbyte - k1->ibyte; ++idx){
            if (k1->csuf[idx] < k2->csuf[idx]){
                return -1;
            }
            else if (k1->csuf[idx] > k2->csuf[idx]){
                return 1;
            }
        }
    }
    return 0;
}

void print_dat(FILE* outf, char* kmer){
    for (int i = 0; i < strlen(kmer); ++i){
        switch(kmer[i]){
            case 'a':
            case 'A':
                fprintf(outf, "A");
            break;
            case 'c':
            case 'C':
                fprintf(outf, "C");
            break;
            case 'g':
            case 'G':
                fprintf(outf, "G");
            break;
            case 't':
            case 'T':
                fprintf(outf, "T");
            break;
            default:
                fprintf(outf, "N");
            break;
        }
    }
    fprintf(outf, "\n");
}

int main(int argc, char* argv[]){
    
    // Set defaults

    char kmers[50][100];
    int kmers_idx = 0;
    char names[50][100];
    int names_idx = 0;
    char readsfile[100];
    char proffile[100];
    char output_prefix[100];   
    readsfile[0] = '\0';
    output_prefix[0] = '\0';
    
    int option_index = 0;
    int ch;
    if (argc == 1){
        help(0);
    }
    
    // Parse arguments
    while ((ch = getopt(argc, argv, "k:n:r:o:h")) != -1){
        switch(ch){
            case 'k':
                strcpy(&kmers[kmers_idx][0], optarg);
                for (int i = strlen(optarg); i < 100; ++i){
                    kmers[kmers_idx][i] = '\0';
                }
                kmers_idx++;
                break;
            case 'h':
                help(0);
                break;
            case 0:
                break;
            case 'o':
                strcpy(&output_prefix[0], optarg);
                break;
            default:
                help(0);
                break;
        }
    }
    // Validate args
    if (strlen(output_prefix) == 0){
        fprintf(stderr, "ERROR: --output_prefix / -o required.\n");
        exit(1);
    }
    if (kmers_idx == 0){
        fprintf(stderr, "ERROR: --kmers / -k is required\n");
        exit(1);
    }
    if (kmers_idx == 1){
        fprintf(stderr, "ERROR: cannot demultiplex with only one species\n");
        exit(1);
    }
    
    // Load k-mer tables
    int num_tables = kmers_idx;
    Kmer_Stream* tables[num_tables];
    // Ensure all k-mers are same
    int kmer = -1;
    for (int i = 0; i < num_tables; ++i){
        fprintf(stderr, "loading %s\n", kmers[i]);
        tables[i] = Open_Kmer_Stream(kmers[i]);
        if (tables[i] == NULL){
            fprintf(stderr, "ERROR opening %s\n", kmers[i]);
            exit(1);
        }
        else{
            if (kmer != -1 && tables[i]->kmer != kmer){
                fprintf(stderr, "ERROR: conflicting k-mer lengths: %d %d\n", kmer, tables[i]->kmer);
                exit(1);
            }
            kmer = tables[i]->kmer;
            fprintf(stderr, "Loaded %d-mer table %s\n", tables[i]->kmer, kmers[i]);
            
        }
    }
    // prepare output files
    FILE* outs[num_tables];
    char namebuf[500];
    for (int i = 0; i < num_tables; ++i){
        sprintf(&namebuf[0], "%s.%d.kmers", output_prefix, i);
        outs[i] = fopen(namebuf, "w");
        if (outs[i] == NULL){
            fprintf(stderr, "ERROR opening %s for writing.\n", &namebuf[0]);
            exit(1);
        }
    }
    for (int i = 0; i < num_tables; ++i){
        First_Kmer_Entry(tables[i]);
    }
    char kmer_text[150]; 
    int ntie = 0;
    int ties[num_tables];
    int finished = 0;
    while (!finished){
        // Each iteration: check for unique k-mers
        // Then increment only iterators with lowest sort-order k-mers
        // A k-mer is unique if it's lower sort order than the other streams
        // and does not match others
        int nvalid = 0;
        for (int i = 0; i < num_tables; ++i){
            if (tables[i]->csuf != NULL){
                nvalid++;
            }
        }
        if (nvalid < 2){
            // Not enough k-mers to compare.
            finished = 1;
            break;
        } 
        else{
            ntie = 0;
            int min_idx = 0;
            while (tables[min_idx]->csuf == NULL){
                ++min_idx;
            }
            // Compare all k-mers
            for (int i = min_idx + 1; i < num_tables; ++i){
                if (tables[i]->csuf != NULL){ 
                    int comp = kcomp(tables[i], tables[min_idx]);
                    if (comp < 0){
                        min_idx = i;
                        ntie = 0;
                    }
                    else if (comp == 0){
                        // tie
                        ties[ntie] = i;
                        ntie++;
                    }
                }
            }
            // If no tie, print the lowest.
            if (ntie == 0){
                char* b = Current_Kmer(tables[min_idx], &kmer_text[0]);
                int c = Current_Count(tables[min_idx]);
                if (c == 1){
                    print_dat(outs[min_idx], b);
                }
            }
            // Only increment lowest-value iterators.
            int incremented = 0;
            if (tables[min_idx]->csuf != NULL){
                Next_Kmer_Entry(tables[min_idx]);
                incremented = 1;
            }
            for (int i = 0; i < ntie; ++i){
                if (tables[ties[i]]->csuf != NULL){
                    Next_Kmer_Entry(tables[ties[i]]);
                    incremented = 1;
                }   
            }
            if (!incremented){
                // Break out of the loop
                finished = 1;
                break;
            }
        }
    }
    
    // Any remaining entries in the last iterator are unique.
    for (int i = 0; i < num_tables; ++i){
        while (tables[i]->csuf != NULL){
            char* b = Current_Kmer(tables[i], &kmer_text[0]);
            print_dat(outs[i], b);
            Next_Kmer_Entry(tables[i]);
        }
    }
    for (int i = 0; i < num_tables; ++i){
        fclose(outs[i]);
    }
    //free(entries);
    //free(bases);
    return 0;

}
