//
//  classifier.c
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#include "classifier.h"
#include "File.h"
#include "getline.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct s_Classifier
{
    LZ78 *models;
    unsigned int nmodels;
    unsigned int models_size;
};

////////////////////////////////////////////////////////////////////////////////
Classifier ClassifierCreate(void)
{
    Classifier classifier = (Classifier)calloc(1, sizeof(struct s_Classifier));
    if (!classifier)
    {
        fprintf(stderr, "%s, %u: memory allocation of %lu bytes failed\n\n", __FILE__, __LINE__, sizeof(struct s_Classifier));
        exit(-1);
    }
    return classifier;
}

////////////////////////////////////////////////////////////////////////////////
void ClassifierDestroy(Classifier classifier)
{
    for(int i=0; i<classifier->nmodels; i++)
        LZ78Destroy(classifier->models[i]);

    free(classifier->models);
    free(classifier);
}

////////////////////////////////////////////////////////////////////////////////
void ClassifierAdd(Classifier classifier, const char *name, const char *file_path, unsigned int max_depth)
{
    LZ78 model = NULL;
    for (unsigned int i = 0; i < classifier->nmodels; i++)
    {
        if (!strcmp(LZ78Name(classifier->models[i]), name))
        {
            model = classifier->models[i];
        }
    }

    if (!model)
    {
        if (classifier->nmodels == classifier->models_size)
        {
            classifier->models_size = (classifier->models_size > 0) ? classifier->models_size * 2 : 2048;
            classifier->models = (LZ78*)realloc(classifier->models, classifier->models_size * sizeof(LZ78));
            if (!classifier->models)
            {
                fprintf(stderr, "%s, %u: memory allocation of %lu bytes failed\n\n", __FILE__, __LINE__, sizeof(LZ78 *) * classifier->models_size);
                exit(-1);
            }
        }
        model = LZ78Create(name, max_depth);
        if(!model) {
            fprintf(stderr, "%s, %d: this should not happen, the model is undefined\n\n", __FILE__, __LINE__);
            exit(-1);
        }
        classifier->models[classifier->nmodels++] = model;
    }

    LZ78Build(model, file_path);
}

////////////////////////////////////////////////////////////////////////////////
void ClassifierBatchAdd(Classifier classifier, const char *name2file, unsigned int max_depth)
{
    FILE *fin = fopen(name2file, "r");
    char *buf = NULL;
    size_t buf_size = 0;

    if (!fin)
    {
        fprintf(stderr, "%s, %u: failed to open file %s for reading\n\n", __FILE__, __LINE__, name2file);
        exit(-1);
    }
    while (getline(&buf, &buf_size, fin) > 0)
    {
        if (strrchr(buf, '\n'))
        {
            *strrchr(buf, '\n') = 0;
        }

        const char *name = buf;
        char *file_path = strchr(buf, '\t');
        if (!file_path)
        {
            fprintf(stderr, "%s, %u: Illegal line in  file %s:\n%s\n\n", __FILE__, __LINE__, name2file, buf);
            exit(-1);
        }
        *(file_path++) = 0;
        ClassifierAdd(classifier, name, file_path, max_depth);
    }
    free(buf);
    fclose(fin);
}

////////////////////////////////////////////////////////////////////////////////
const char *ClassifierPredict(Classifier classifier, const char *name, const char *file_path, FILE *fout)
{
    unsigned int i, best_i = 0;
    double best_score = 10000;
    if (classifier->nmodels == 0)
        return NULL;

    File file = FileCreate(file_path);

    fprintf(fout, "%s", name);
    for (i = 0; i < classifier->nmodels; i++)
    {
        double score = LZ78AverageLogScore(classifier->models[i], file);
        fprintf(fout, "\t%lf", score);
        if(score < best_score) {
            best_score = score;
            best_i = i;
        }
    }
    fprintf(fout, "\t%s\n", LZ78Name(classifier->models[best_i]));
    FileDestroy(file);

    return LZ78Name(classifier->models[best_i]);
}

////////////////////////////////////////////////////////////////////////////////
void ClassifierPrintHeader(Classifier classifier, FILE *fout)
{
    fprintf(fout, "Genome_name");

    for (unsigned int i = 0; i < classifier->nmodels; i++)
    {
        fprintf(fout, "\t%s", LZ78Name(classifier->models[i]));
    }
    fprintf(fout, "\tBest_hit\n");
}

////////////////////////////////////////////////////////////////////////////////
void ClassifierPrintStats(Classifier classifier, FILE *fout)
{
    fprintf(fout, "\nNumber of models: %u\n", classifier->nmodels);
    fprintf(fout, "some stats for each model:\n\n");

    for (unsigned int i = 0; i < classifier->nmodels; i++)
    {
        fprintf(fout, "--------------------------------------------------\n");
        LZ78WriteStats(classifier->models[i], fout);
    }
    fprintf(fout, "\n");
}
