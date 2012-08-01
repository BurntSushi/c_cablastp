#ifndef __CABLASTP_ALIGN_H__
#define __CABLASTP_ALIGN_H__

/*
This is adapted from Dan Kortschak's Needleman-Wunsch algorithm from the biogo
package: code.google.com/p/biogo.

It's mostly copied from its original form, but it is optimized specifically
for cablastp to limit allocations and to absolve the need for the biogo/seq.Seq
type.
*/

#include <stdint.h>

#define CABLASTP_ALIGN_SEQ_SIZE 10000

int32_t
cbp_align_ungapped(int32_t window_size, int32_t kmer_size, int32_t id_threshold,
                   char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend);

int32_t
cbp_align_identity(char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend);

struct cbp_align_nw_memory {
    int32_t *table;
    int32_t *zeroes;
    char *ref;
    char *org;
};

struct cbp_align_nw_memory *
cbp_align_nw_memory_init();

void
cbp_align_nw_memory_free(struct cbp_align_nw_memory *mem);

struct cbp_alignment {
    char *ref;
    char *org;
    int32_t length;
};

struct cbp_alignment
cbp_align_nw(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend,
             char *oseq, int32_t ostart, int32_t oend);

int32_t
cbp_align_length_nogaps(char *residues);

#endif
