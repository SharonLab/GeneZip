//
//  File.h
//  
//  A simple file reader. Reads the contents of the file into a buffer, then returns the contents
//  line by line. Initial buffer size is FILE_BUFFER_SIZE, should fit most prokaryotic genome sizes.
//  Buffer size increases as required
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#ifndef File_h
#define File_h
#include <stdio.h>

#define	FILE_BUFFER_SIZE	10000000

typedef struct s_File* File;

// Will load the file into a buffer
File FileCreate(const char* fname);

// Delete file buffer and the object
void FileDestroy(File file);

// Read the next line into the buffer pointed by buf. Will change buf's size if necessary
// and store the new size in bufsize. Returns the number of characters read
size_t FileGetline(char** buf, size_t* bufsize, File file);

// Go back to the beginning of the file
void FileRollBack(File file);

#endif /* File_h */
