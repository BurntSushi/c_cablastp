#ifndef __CABLASTP_COARSE_H__
#define __CABLASTP_COARSE_H__

/* Apparently this is required to make pthread_rwlock* stuff available. */
#define __USE_UNIX98

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>

#include "ds.h"

/* #include "seeds.h" */
#include "seq.h"

struct cbp_link_to_compressed {
    int32_t org_seq_id;
    int16_t coarse_start;
    int16_t coarse_end;
    struct cbp_link_to_compressed *next;
};

struct cbp_link_to_compressed *
cbp_link_to_compressed_init(int32_t org_seq_id,
                            int16_t coarse_start, int16_t coarse_end);

void
cbp_link_to_compressed_free(struct cbp_link_to_compressed *link);

struct cbp_coarse_seq {
    int32_t id;
    struct cbp_seq *seq;
    struct cbp_link_to_compressed *links;
    pthread_rwlock_t lock_links;
};

struct cbp_coarse_seq *
cbp_coarse_seq_init(int32_t id, char *residues, int32_t start, int32_t end);

void
cbp_coarse_seq_free(struct cbp_coarse_seq *seq);

void
cbp_coarse_seq_addlink(struct cbp_coarse_seq *seq,
                       struct cbp_link_to_compressed *newlink);

struct cbp_coarse {
    struct DSVector *seqs;
    /* struct cbp_seeds *seeds; */
    FILE *file_fasta;
    FILE *file_seeds;
    FILE *file_links;
    pthread_rwlock_t lock_seq;
};

struct cbp_coarse *
cbp_coarse_init(int32_t seed_size,
                FILE *file_fasta, FILE *file_seeds, FILE *file_links);

void
cbp_coarse_free(struct cbp_coarse *coarse_db);

struct cbp_coarse_seq *
cbp_coarse_add(struct cbp_coarse *coarse_db,
               char *residues, int32_t start, int32_t end);

struct cbp_coarse_seq *
cbp_coarse_get(struct cbp_coarse *coarse_db, int32_t i);

void
cbp_coarse_save_plain(struct cbp_coarse *coarse_db);

#endif
