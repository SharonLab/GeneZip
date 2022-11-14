//
// cli.h
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#ifndef CLI_H
#define CLI_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lz78.h"

#define DEFAULT_MAX_DEPTH 13

typedef struct s_Usage* Usage;

Usage get_usage(int argc, const char *argv[], const char*);
const char* get_training_name2file_file(Usage);
const char* get_prediction_name2file_file(Usage);
const char* get_out_file(Usage);
unsigned int get_max_depth(Usage);
void destroy_usage(Usage);

#endif //CLI_H
