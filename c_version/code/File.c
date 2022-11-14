//
//  File.c
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#include "File.h"
#include <stdlib.h>
#include <string.h>

struct s_File {
    char* buf;
    size_t bufsize;
    size_t curr;
};

////////////////////////////////////////////////////////////////////////////////
File FileCreate(const char* fname) {
    FILE* fin = fopen(fname, "rb");
    if(!fin) {
        fprintf(stderr, "%s, %d: failed to read %s\n\n", __FILE__, __LINE__, fname);
        exit(-1);
    }

    File file = (File)malloc(sizeof(struct s_File));
    if(!file) {
        fprintf(stderr, "%s, %d: failed allocate %lu bytes\n\n", __FILE__, __LINE__, sizeof(struct s_File));
        exit(-1);
    }

    file->bufsize = FILE_BUFFER_SIZE;
    file->buf = (char*)calloc(file->bufsize, sizeof(char));
    if(!(file->buf)) {
        fprintf(stderr, "%s, %d: failed allocate %lu bytes\n\n", __FILE__, __LINE__, sizeof(char)*file->bufsize);
        exit(-1);
    }
    file->curr = 0;
    size_t n = 0, nread = 0;

    while((nread = fread(file->buf+n, 1, FILE_BUFFER_SIZE, fin)) == FILE_BUFFER_SIZE)  {
        file->bufsize += FILE_BUFFER_SIZE;
        file->buf = (char*)realloc(file->buf, file->bufsize*sizeof(char));
        if(!(file->buf)) {
            fprintf(stderr, "%s, %d: failed allocate %lu bytes\n\n", __FILE__, __LINE__, sizeof(char)*file->bufsize);
            exit(-1);
        }
        n += nread;
    }
    n += nread;
    // We must have at least one byte left in the buffer because of the condition we put in the while loop
    file->buf[n] = 0;
    if(feof(fin)) {
        fclose(fin);
    } else if(ferror(fin)) {
        if(fin) { fclose(fin); }
        fprintf(stderr, "%s, %d: file reading ended due to an error for file %s\n\n", __FILE__, __LINE__, fname);
        exit(-1);
    }

    return file;
}

////////////////////////////////////////////////////////////////////////////////
void FileDestroy(File file) {
    free(file->buf);
    free(file);
}

////////////////////////////////////////////////////////////////////////////////
void FileRollBack(File file) {
    file->curr = 0;
}

////////////////////////////////////////////////////////////////////////////////
size_t FileGetline(char** buf, size_t* bufsize, File file) {
    // Skip empty lines
    while(file->buf[file->curr] == '\n')
        file->curr++;

    if(file->buf[file->curr] == 0)
        return 0;

    size_t n = file->curr;
    while(file->buf[n] && (file->buf[n] != '\n'))
        n++;

    // If the line ends with \n then we need to include it in the returned line. Else, we are at the
    // end of the file
    size_t line_len = n - file->curr + (file->buf[n] == '\n');
    if(line_len+1 > *bufsize) {
        *bufsize = line_len+1;
        *buf = (char*)realloc(*buf, (*bufsize)*sizeof(char));
        if(!(*buf)) {
            fprintf(stderr, "%s, %d: failed allocate %lu bytes\n\n", __FILE__, __LINE__, (*bufsize)*sizeof(char));
            exit(-1);
        }
    }
    strncpy(*buf, file->buf+file->curr, line_len);
    (*buf)[line_len] = 0;

    // file->buf[n] is either \n or \0. If it is \n then we want file->curr to point to the next char.
    // else we want it to be n (reached the end of the file)
    file->curr = (file->buf[n] == '\n')? n+1 : n;

    return line_len;
}
