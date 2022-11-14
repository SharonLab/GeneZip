//
//  lz78.h
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//
#ifndef lz78_h
#define lz78_h

#define	MAX_DEPTH	17

#include "File.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

typedef struct s_LZ78 *LZ78;

// Constructor: name is the name of the model, max_depth is the maximum
// allowed depth of any path in the tree, root to leaf 
LZ78	LZ78Create(const char *name, unsigned int max_depth);

// Destructor
void	LZ78Destroy(LZ78 lz);

// Build the actual model
size_t	LZ78Build(LZ78 lz, const char *file);

// Calculate the average log-loss for the sequences in file (considered together)
double	LZ78AverageLogScore(LZ78 lz, File file);

// Get the name of the model
const char*	LZ78Name(LZ78 lz);

// Returns the number of inner nodes in the tree
unsigned int	LZ78NumInnerNodes(LZ78 lz);

// Returns the maximum depth with all inner nodes present
unsigned int	LZ78MaxCompleteDepth(LZ78 lz);

// Returns the longest path from root to leaf in the tree
unsigned int	LZ78LongestPathRootToLeaf(LZ78 lz);

// Writes statistics about the model
void	LZ78WriteStats(LZ78 lz, FILE* fout);

#endif /* lz78_h */
