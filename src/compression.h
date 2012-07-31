#ifndef __CABLASTP_COMPRESSION_H__
#define __CABLASTP_COMPRESSION_H__

/* Apparently this is required to make pthread_rwlock* stuff available. */
#define __USE_UNIX98

#include <pthread.h>
#include <stdint.h>

#include "ds.h"

#include "align.h"
#include "coarse.h"
#include "compressed.h"
#include "database.h"
#include "seq.h"

struct cbp_compress_workers {
    pthread_t *threads;
    int32_t num_workers;
    struct DSQueue *jobs;
    void *args;
};

struct cbp_compress_workers *
cbp_compress_start_workers(struct cbp_database *db, int32_t num_workers);

void
cbp_compress_join_workers(struct cbp_compress_workers *workers);

void
cbp_compress_free_workers(struct cbp_compress_workers *workers);

void
cbp_compress_send_job(struct cbp_compress_workers *workers,
                      struct cbp_seq *org_seq);

struct cbp_compressed_seq *
cbp_compress(struct cbp_coarse *coarse_db, struct cbp_seq *org_seq,
             struct cbp_align_nw_memory *mem);

#endif
