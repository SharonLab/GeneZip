#ifndef GETLINE_H
#define GETLINE_H
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

ssize_t getdelim(char** buf, size_t* bufsiz, int delimiter, FILE* fp);
ssize_t getline(char** buf, size_t* bufsiz, FILE* fp);

#endif
