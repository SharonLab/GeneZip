//
//  classifier.h
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#ifndef classifier_h
#define classifier_h

#include "lz78.h"

typedef struct s_Classifier* Classifier;

Classifier ClassifierCreate(void);

void ClassifierDestroy(Classifier classifier);

// Add one file to one model (or create the model if it does not exist yet)
void ClassifierAdd(Classifier classifier, const char* name, const char* file_path, unsigned int max_depth);

// Add multiple files at once to build the model. name2file contains, in each line, the name of the model and a path to the  file
void ClassifierBatchAdd(Classifier classifier, const char* name2file_path, unsigned int max_depth);

// Does the predictions, writes the values of each model to an output file (fout)
const char* ClassifierPredict(Classifier classifier, const char* name, const char* file_path, FILE* fout);

// Prints the  header of the output file which contains the order of the models
void ClassifierPrintHeader(Classifier classifier, FILE* fout);

// Prints stats about the different models
void ClassifierPrintStats(Classifier classifier, FILE* fout);

#endif /* classifier_h */

