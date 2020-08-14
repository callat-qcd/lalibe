#ifndef __recursive_coloring_h__
#define __recursive_coloring_h__

//This is Andrea's file that has been modified slightly for chroma.
//I put many declarations in a header file.
//Arjun Singh Gambhir

//I add this here so I can cal some Layout stuff later on.
namespace Chroma
{

//I put the struct typedef in the header.
/************************************************************************/
typedef struct meshVars {
   unsigned int d;		/* Number of dimensions */
   unsigned int *ptsPerDim;     /* points per dimension */
   unsigned int dmax;		/* which dimension has most points */
   unsigned int twoTod;         /* 2^d */
   unsigned int *logPts;        /* log2(ptsPerDim) */
   char *RBorder;  /* Red-black order of a d-dim torus with 2 pts per dim */
   char *PTbits;   /* Bits for d coordinates of a point (d x maxlog2(Pts)) */
} meshVars;

/************************************************************************/

/* Turn integer to binary (chars 0 or 1). log(n) should be <= arraySize */
void int2bin(unsigned int n, char *array, unsigned int arraySize) ;

/* Convert binary string to 10-base int. Assumes array has at least arraySize */
unsigned int bin2int(char *array, unsigned int arraySize);

/************************************************************************/
/* RBorder array: Permutation array of 2^d nodes.
 * However, permutation nodes are kept in binary (char)
 * so it is a one d array which can be viewed as:
 *         d columns
 *  2^d   <11000110>   -> represents a number from 0 to 2^d-1
 *  rows  <01001010>  
 *        <10111010>  
 */

void findRBorder(struct meshVars *mesh);

/************************************************************************/

void hierOrderSetUp(struct meshVars *mesh);

/************************************************************************/

unsigned int hierOrderPoint(unsigned int *coord, struct meshVars *mesh);

/************************************************************************/
/* Deallocate mesh info */
void freeMeshVars(struct meshVars *mesh);
/************************************************************************
 * Hadamard vector related functions 
 ************************************************************************/

/* hadaColPerm produces the index sequence of hadamard columns which 
 * correspond to a red-black ordering of a power of two dimension 
 *
 * **CAUTION** IT ONLY WORKS FOR N power of 2 **
 *
 * Input: 
 * 	unsigned int N: dimension the vector space 
 * 	unsigned int Hpsize: number of indices needed in the sequence
 * Output: 
 * 	unsigned int *Hperm: the array containing the indices
 */
void hadaColPerm(unsigned int N, unsigned int* Hperm, unsigned int Hpsize);
/************************************************************************/
unsigned int num1bits(unsigned int a);

int Hada_element(unsigned int i, unsigned int j);
/************************************************************************
 * PROBLEM SPECIFIC FUNCTIONS
 ************************************************************************/

/************************************************************************/
void index2coord(unsigned int i, struct meshVars *mesh, unsigned int *coord);
/************************************************************************/
/* Create permutation array for all i rows local on this processor */
/* This is problem specific and SHOULD BE ADAPTED */
void hierPerm(struct meshVars *mesh, unsigned int *perm, unsigned int N);

} //Chroma namespace

#endif
