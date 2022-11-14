//
// cli.c
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#include "cli.h"

struct s_Usage
{
    unsigned int max_depth;
    const char* training_name2file_file;
    const char* prediction_name2file_file;
    const char* out_file;
};

static void print_help(const char* prog, const char* version) {
    fprintf(stderr, "GeneZip, %s\n", version);
    fprintf(stderr, "\nUsage: %s -i <training-name2file> -t <predict-name2file> -o <output> [-d <max-depth>]\n\n", prog);
    fprintf(stderr, "  <training-name2file>: a file with the list of fasta files for the cluster models in the format\n");
    fprintf(stderr, "                        <cluster-name>\t<fasta-file>\n");
    fprintf(stderr, "  <predict-name2file> : a file with the list of fasta files for prediction, format is\n");
    fprintf(stderr, "                        <cluster-name>\t<fasta-file>\n");
    fprintf(stderr, "  <max-depth>         : maximum depth allowed for the context tree, between 1 and %d (default: %u)\n", MAX_DEPTH, DEFAULT_MAX_DEPTH);
    fprintf(stderr, "  <out-file>          : name of the output file\n\n");
}

Usage get_usage(int argc, const char *argv[], const char* version) {
    Usage usage = (Usage) malloc(sizeof(struct s_Usage));
    if(!usage) {
        fprintf(stderr, "%s, %d: failed allocate %lu bytes\n\n", __FILE__, __LINE__, sizeof(struct s_Usage));
        exit(-1);
    }

    usage->max_depth = DEFAULT_MAX_DEPTH;
    usage->training_name2file_file = NULL;
    usage->prediction_name2file_file = NULL;
    usage->out_file = NULL;

    if ((argc == 1) || ((argc == 2) && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || !strcmp(argv[1], "--help")))) {
        print_help(argv[0], version);
        exit(0);
    }
    for(int i=1; i<argc-1; i+=2) {
        if(!strcmp(argv[i], "-i")) {
            usage->training_name2file_file  = argv[i+1];
            FILE* fp = fopen(usage->training_name2file_file, "r");
            if(!fp) {
                fprintf(stderr, "Error (-i): failed to open training file (%s) for reading\n\n", usage->training_name2file_file );
                exit(1);
            }
            fclose(fp);
        }
        else if(!strcmp(argv[i], "-t")) {
            usage->prediction_name2file_file = argv[i+1];
            FILE* fp = fopen(usage->prediction_name2file_file, "r");
            if(!fp) {
                fprintf(stderr, "Error (-t): failed to open prediction file (%s) for reading\n\n", usage->prediction_name2file_file);
                exit(1);
            }
            fclose(fp);
        }
        else if(!strcmp(argv[i], "-o")) {
            usage->out_file = argv[i+1];
        }
        else if(!strcmp(argv[i], "-d")) {
            usage->max_depth = (unsigned int)atoi(argv[i+1]);
            if ((usage->max_depth == 0) || (usage->max_depth > MAX_DEPTH)) {
                fprintf(stderr, "Error (-d): illegal max_depth (%s), must be an integer between 1 to %d\n\n", argv[i+1], MAX_DEPTH);
                exit(1);
            }
        }
        else {
            fprintf(stderr, "Error: unknown option %s\n", argv[i]);
            print_help(argv[0], version);
            exit(1);
        }
    }

    if(!usage->training_name2file_file || !usage->prediction_name2file_file || !usage->out_file) {
        fprintf(stderr, "\nError: training (-i), prediction (-t) or output (-o) files are not specified\n");
        print_help(argv[0], version);
        exit(1);
    }

    return usage;
}

const char* get_training_name2file_file(Usage usage) {
    return usage->training_name2file_file;
}

const char* get_prediction_name2file_file(Usage usage) {
    return usage->prediction_name2file_file;
}

const char* get_out_file(Usage usage) {
    return usage->out_file;
}

unsigned int get_max_depth(Usage usage) {
    return usage->max_depth;
}

void destroy_usage(Usage usage) {
    free(usage);
}
