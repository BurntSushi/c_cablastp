#ifndef __CABLASTP_COMPRESSED_H__
#define __CABLASTP_COMPRESSED_H__

#include <stdint.h>
#include <stdio.h>

#include "ds.h"

struct cbp_link_to_coarse {
    char *diff;
    int32_t coarse_seq_id;
    int16_t coarse_start;
    int16_t coarse_end;
    struct cbp_link_to_coarse *next;
};

struct cbp_link_to_coarse *
cbp_link_to_coarse_init(int32_t coarse_seq_id,
                        int16_t coarse_start, int16_t coarse_end,
                        char *align_ref, char *align_tgt);

struct cbp_link_to_coarse *
cbp_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               int16_t coarse_start, int16_t coarse_end);

void
cbp_link_to_coarse_free(struct cbp_link_to_coarse *link);

struct cbp_compressed_seq {
    int32_t id;
    char *name;
    struct cbp_link_to_coarse *links;
};

struct cbp_compressed_seq *
cbp_compressed_seq_init(int32_t id, char *name);

void
cbp_compressed_seq_free(struct cbp_compressed_seq *seq);

void
cbp_compressed_seq_addlink(struct cbp_compressed_seq *seq,
                           struct cbp_link_to_coarse *link);

struct cbp_compressed {
    struct DSVector *seqs;
    FILE *file_compressed;
    FILE *file_index;
};

struct cbp_compressed *
cbp_compressed_init(FILE *file_compressed, FILE *file_index);

void
cbp_compressed_free(struct cbp_compressed *com_db);

int32_t
cbp_compressed_size(struct cbp_compressed *com_db);

void
cbp_compressed_add(struct cbp_compressed *com_db,
                   struct cbp_compressed_seq *seq);

void
cbp_compressed_save_plain(struct cbp_compressed *com_db);

#endif
