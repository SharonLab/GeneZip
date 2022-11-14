//
//  lz78.c
//
//  Created by Or Leibovich, Yochai Meir, and Itai Sharon, last updated on 11/Sep/22
//

#include "lz78.h"
#include <assert.h>

#define	CHECK_BIT(idx)	((lz->mem[(idx) >> 3] & (128 >> ((idx) & 7))) > 0)

////////////////////////////////////////////////////////////////////////////////
// LZ78
////////////////////////////////////////////////////////////////////////////////
struct s_LZ78
{
    unsigned char *mem;		// Bit array that keeps the inner nodes
    unsigned int mem_size;	// Memory size allocated for mem
    unsigned int max_depth;	// Maximum depth including leaves, specified by the user. 
    unsigned int leaf_count;	// Total number of leaves (paths) in the tree. Used for calculating log-loss
    unsigned int full_depth;	// Maximum depth in which all inner nodes are present
    char* name;
    size_t num_nodes_in_depth[MAX_DEPTH+1];	// Keeps the number of inner nodes in each depth
};

////////////////////////////////////////////////////////////////////////////////
void    LZ78WriteStats(LZ78 lz, FILE* fout) {
    fprintf(fout, "Name:                      %s\n", LZ78Name(lz));
    fprintf(fout, "Node array size:           %u\n", lz->mem_size);
    fprintf(fout, "Number of inner nodes:     %u\n", LZ78NumInnerNodes(lz));
    fprintf(fout, "Max complete depth:        %u\n", LZ78MaxCompleteDepth(lz));
    fprintf(fout, "Longest path (root->leaf): %u\n", LZ78LongestPathRootToLeaf(lz));
    fprintf(fout, "Number of inner node in each depth (%c of possible nodes):\n", '%');
    fprintf(fout, "Depth\tNNodes\tNFull\t%c of full\n", '%');

    fprintf(fout, "0\t1\t1\t100.0\n");
    for(int i=1; i<lz->max_depth; i++) {
       size_t full_n = ((size_t)4) << 2*(i-1);
       fprintf(fout, "%d\t%lu\t%lu\t%.1lf\n", i, lz->num_nodes_in_depth[i], full_n, 100.0*((double)lz->num_nodes_in_depth[i])/((double) full_n));
    }
    fprintf(fout, "\nNumber of leaves:\t%u\n", lz->leaf_count);
}

////////////////////////////////////////////////////////////////////////////////
// len_bases contains the number of bits required for each depth. At each depth, we need 4 nodes for each leaf
// in the previous depth. depth 1: 4, depth 2: 4+16, depth 3: 4+16+64, etc. Maximum depth is 17
size_t len_bases[MAX_DEPTH+1] = {}; // {0, 4, 20, 84, 340, 1364, 5460, 21844, 87380, 349524, 1398100, 5592404, 22369620, 89478484, 357913940, 1431655764, 5726623060, 22906492244};

// Returns the memory size required for keeping all inner nodes up to depth max_depth,
// where max_depth may contain only leaves
size_t calc_mem_size(unsigned int max_depth)
{
    return (len_bases[max_depth-1] / 8) + 1;
}

////////////////////////////////////////////////////////////////////////////////
void AddNode(LZ78 lz, size_t node_index, unsigned int depth)
{
    lz->mem[node_index >> 3] |= (128 >> (node_index & 7));
    lz->num_nodes_in_depth[depth]++;

    // One leaf became an inner node, 4 new leaves were added, 3 new leaves added in total
    lz->leaf_count += 3;
}

////////////////////////////////////////////////////////////////////////////////
LZ78 LZ78Create(const char *name, unsigned int depth)
{
    if((depth == 0) || (depth > MAX_DEPTH)) {
        fprintf(stderr, "%s, %d: requested depth (%u) is illegal, must be between 1 to %d. Aborting\n\n", __FILE__, __LINE__, depth, MAX_DEPTH);
        exit(-1);
    }

    LZ78 lz = (LZ78)calloc(1, sizeof(struct s_LZ78));
    if (!lz)
    {
        fprintf(stderr, "%s, %d: memory allocation of %lu bytes failed\n\n", __FILE__, __LINE__, sizeof(struct s_LZ78));
        exit(-1);
    }

    lz->name = strdup(name);
    if(!lz->name) {
        fprintf(stderr, "%s, %d: str duplication failed\n\n", __FILE__, __LINE__);
        exit(-1);
    }

    // len_bases is the starting index for each depth i AND the number of indices
    // required for depth i-1     
    for(int i=1; i<=MAX_DEPTH; i++)
        len_bases[i] = len_bases[i-1] + (4 << 2*(i-1));

    lz->mem_size = calc_mem_size(depth);
    lz->max_depth = depth;
    lz->mem = (unsigned char *)calloc(lz->mem_size, sizeof(unsigned char));
    if (!lz->mem)
    {
        fprintf(stderr, "%s, %d: memory allocation of %lu bytes failed\n\n", __FILE__, __LINE__, lz->mem_size*sizeof(unsigned char));
        exit(-1);
    }

    lz->num_nodes_in_depth[0] = 1;
    lz->leaf_count = 4;
    lz->full_depth = 0;

    return lz;
}

////////////////////////////////////////////////////////////////////////////////
void LZ78Destroy(LZ78 lz) {
    free(lz->mem);
    free(lz->name);
    free(lz);
}

////////////////////////////////////////////////////////////////////////////////
// For a sequence s=s1s2...sn, the index in the bit memory can be calculated 
// using m(s1..si) = (4^0+..+4^(i-1))-1 + 4*m(s1..si-1) 
// There is some bug here
size_t LZ78Build(LZ78 lz, const char *file_path)
{
    File file = FileCreate(file_path);

    char *buf = NULL;
    size_t buf_size = 0;
    size_t n;

    unsigned int curr_depth = 1; // len is at least 1
    size_t curr_sequence = 0;  // sequence value.

    while ((n = FileGetline(&buf, &buf_size, file)) != 0) // EOF)
    {
        if ((*buf == '>') || (*buf == '\n'))
        {
            curr_depth = 1;
            curr_sequence = 0;
            continue;
        }
        char *p = buf, *pe = buf + n;

        // buf may or may not end with \n
	if(*(pe-1) == '\n')
		pe--;

        for (; p < pe; p++)
        {
            // toupper
            if (*p >= 'a')
                (*p) -= 32;

            // This will start a new path
            if (*p == 'N')
            {
                curr_depth = 1;
                curr_sequence = 0;
                continue;
            }

            // For (*p >> 1) this is what we get:
            // A = 0b100000
            // C = 0b100001
            // G = 0b100011
            // T = 0b101010
            // Order will be ACTG
            int i = (*p >> 1) & 3;
            curr_sequence |= i;

            // As far as I understand this can only happen if lz->max_depth == 1
            if(curr_depth > lz->max_depth-1) {
                curr_depth = 1;
                curr_sequence = 0;
                continue;
            }

            size_t current_index = len_bases[curr_depth - 1] + curr_sequence;
            assert((current_index >> 3) < lz->mem_size);
            if (!CHECK_BIT(current_index))
            {
                AddNode(lz, current_index, curr_depth);
                curr_depth = 1;
                curr_sequence = 0;
                continue;
            }

            if(curr_depth == lz->max_depth-1) {
                curr_depth = 1;
                curr_sequence = 0;
            }
            else
            {
                curr_sequence <<= 2;
                curr_depth++;
            }
        }
    }

    // Check what is the maximum level with all nodes
    lz->full_depth = 0;
    while((lz->full_depth+1 < lz->max_depth) && (lz->num_nodes_in_depth[lz->full_depth+1] == (4 << 2*(lz->full_depth))))
        lz->full_depth++;
    
    free(buf);
    FileDestroy(file);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
double LZ78AverageLogScore(LZ78 lz, File file)
{
    char *buf = NULL;
    size_t buf_size = 0;
    size_t n;

    unsigned int nchars = 0;
    unsigned int actual_nchars = 0;
    unsigned int leaf_count = 0;

    unsigned int curr_depth = 1;
    size_t curr_sequence = 0;

    FileRollBack(file);

    while ((n = FileGetline(&buf, &buf_size, file)) != 0) // EOF)
    {
        if ((*buf == '>') || (*buf == '\n'))
        {
            curr_depth = 1;
            curr_sequence = 0;
            continue;
        }
        char *p = buf, *pe = buf + n;

        // buf may or may not end with \n
	if(*(pe-1) == '\n')
		pe--;

        for (; p < pe; p++)
        {
            if (*p >= 'a')
                (*p) -= 32;
            if (*p == 'N')
            {
                curr_depth = 1;
                curr_sequence = 0;
                continue;
            }
            // For (*p >> 1) this is what we get:
            // A = 0b100000
            // C = 0b100001
            // G = 0b100011
            // T = 0b101010
            // Order will be ACTG
            int i = (*p >> 1) & 3;
            nchars++;
            curr_sequence |= i;
            size_t current_index = len_bases[curr_depth - 1] + curr_sequence;

            // The first part is an optimization: no need to check if the node 
            // exists for depths with all inner nodes present 
            if ((curr_depth <= lz->full_depth) || ((curr_depth < lz->max_depth) && CHECK_BIT(current_index)))
            {
                curr_sequence <<= 2;
                curr_depth++;
            }
            else
            {
                leaf_count++;
                curr_depth = 1;
                curr_sequence = 0;
                actual_nchars = nchars;
                continue;
            }
        }
    }
    free(buf);
    return (log(lz->leaf_count) / log(2) * leaf_count) / actual_nchars; // log_2
}

////////////////////////////////////////////////////////////////////////////////
const char *LZ78Name(LZ78 lz)
{
    return lz->name;
}

////////////////////////////////////////////////////////////////////////////////
unsigned int LZ78NumInnerNodes(LZ78 lz)
{
    size_t n = 0;
    for(int i=0; i<lz->max_depth; i++)
        n += lz->num_nodes_in_depth[i];

    return n;
}
////////////////////////////////////////////////////////////////////////////////
unsigned int LZ78MaxCompleteDepth(LZ78 lz)
{
    return lz->full_depth;
}

////////////////////////////////////////////////////////////////////////////////
unsigned int LZ78LongestPathRootToLeaf(LZ78 lz)
{
    unsigned int i;
    for(i=0; i<lz->max_depth; i++)
        if(lz->num_nodes_in_depth[i] == 0)
            break;
    return i;
}
