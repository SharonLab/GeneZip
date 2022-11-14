//
//  main.c
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "getline.h"
#include "classifier.h"
#include "cli.h"

#define	VERSION	"v1.00, 11/Sep/2022"

/////////////////////////////////////////////////////////////////////////////////////////////////////
int run_lz_classifier(Usage usage) {
    char *buf = NULL;
    size_t buf_size = 0;
    Classifier classifier;
    time_t rawtime;
    struct tm *timeinfo;
    char time_buf[1024];

    // Open the output stream
    FILE* output_stream = fopen(get_out_file(usage), "wt");
    if(!output_stream) {
        fprintf(stderr, "Error: Cannot create output file %s\n\n", get_out_file(usage));
        return 1;
    }

    fprintf(stderr, "GeneZip, %s\n", VERSION);

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(time_buf, "%s", asctime(timeinfo));
    *strchr(time_buf, '\n') = 0;
    fprintf(stderr, "%s\tStarting\n", time_buf);

    classifier = ClassifierCreate();

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(time_buf, "%s", asctime(timeinfo));
    *strchr(time_buf, '\n') = 0;
    fprintf(stderr, "%s\tTraining\n", time_buf);

    ClassifierBatchAdd(classifier, get_training_name2file_file(usage), get_max_depth(usage));

    FILE *fin = fopen(get_prediction_name2file_file(usage), "r");
    if (!fin)
    {
        fprintf(stderr, "\nError: failed to read file %s\n", get_prediction_name2file_file(usage));
        return 1;
    }

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(time_buf, "%s", asctime(timeinfo));
    *strchr(time_buf, '\n') = 0;
    fprintf(stderr, "%s\tPredicting\n", time_buf);

    ClassifierPrintHeader(classifier, output_stream);
    while (getline(&buf, &buf_size, fin) != EOF)
    {
        if (strrchr(buf, '\n'))
        {
            *strrchr(buf, '\n') = 0;
        }
        const char *name = buf;
        char *file_path = strchr(buf, '\t');
        if (!file_path)
        {
            fprintf(stderr, "%s, %u: Illegal line in file %s:\n%s\n\n", __FILE__, __LINE__,  get_prediction_name2file_file(usage), buf);
            exit(-1);
        }
        *(file_path++) = 0;

        ClassifierPredict(classifier, name, file_path, output_stream);
    }
    fclose(fin);
    free(buf);
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(time_buf, "%s", asctime(timeinfo));
    *strchr(time_buf, '\n') = 0;
    fprintf(stderr, "%s\tDone\n", time_buf);
    ClassifierPrintStats(classifier, stderr);
    ClassifierDestroy(classifier);

    fclose(output_stream);

    return 0;
}

int main(int argc, const char *argv[])
{
    Usage usage = get_usage(argc, argv, VERSION);
    int results = run_lz_classifier(usage);
    destroy_usage(usage);

    return results;
}
