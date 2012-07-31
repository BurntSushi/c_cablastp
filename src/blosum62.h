#ifndef __CABLASTP_BLOSUM62_H__
#define __CABLASTP_BLOSUM62_H__

#include <stdint.h>

#define BLOSUM62_ALPHABET "ABCDEFGHIKLMNPQRSTVWXYZ"

#define BLOSUM62_SIZE 24

int32_t
blosum62_residue_to_index(char residue);

const int BLOSUM62_MATRIX[BLOSUM62_SIZE][BLOSUM62_SIZE];

#endif
