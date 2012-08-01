#include <stdio.h>
#include <stdlib.h>

#include "compression.h"
#include "flags.h"

struct worker_args {
    struct cbp_database *db;
    struct DSQueue *jobs;
};

struct extend_match {
    int32_t rlen;
    int32_t olen;
};

struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend,
             char *oseq, int32_t ostart, int32_t oend);

static int32_t
add_without_match(struct cbp_coarse *coarse_db,
                  struct cbp_seq *org_seq, int32_t ostart, int32_t oend);

static void *
cbp_compress_worker(void *data);

static int32_t
min(int32_t a, int32_t b);

struct cbp_compress_workers *
cbp_compress_start_workers(struct cbp_database *db, int32_t num_workers)
{
    struct DSQueue *jobs;
    struct cbp_compress_workers *workers;
    struct worker_args *wargs;
    int32_t i, errno;

    jobs = ds_queue_create(20);

    wargs = malloc(sizeof(*wargs));
    assert(wargs);
    wargs->db = db;
    wargs->jobs = jobs;

    workers = malloc(sizeof(*workers));
    assert(workers);
    workers->threads = malloc(num_workers * sizeof(*workers->threads));
    assert(workers->threads);
    workers->num_workers = num_workers;
    workers->jobs = jobs;
    workers->args = (void*) wargs;

    for (i = 0; i < num_workers; i++) {
        errno = pthread_create(&workers->threads[i], NULL,
            cbp_compress_worker, (void*) wargs);
        if (errno != 0) {
            fprintf(stderr,
                "cbp_compress_start_workers: Could not start thread. Errno: %d",
                errno);
            exit(1);
        }
    }

    return workers;
}

void
cbp_compress_join_workers(struct cbp_compress_workers *workers)
{
    int i, errno;

    ds_queue_close(workers->jobs);

    for (i = 0; i < workers->num_workers; i++) {
        errno = pthread_join(workers->threads[i], NULL);
        if (errno != 0) {
            fprintf(stderr,
                "cbp_compress_join_workers: Could not join thread. Errno: %d",
                errno);
            exit(1);
        }
    }
}

void
cbp_compress_free_workers(struct cbp_compress_workers *workers)
{
    ds_queue_free(workers->jobs);
    free(workers->args);
    free(workers->threads);
    free(workers);
}

void
cbp_compress_send_job(struct cbp_compress_workers *workers,
                      struct cbp_seq *org_seq)
{
    ds_queue_push(workers->jobs, (void*) org_seq);
}

struct cbp_compressed_seq *
cbp_compress(struct cbp_coarse *coarse_db, struct cbp_seq *org_seq,
             struct cbp_align_nw_memory *mem)
{
    struct extend_match mlens;
    struct cbp_coarse_seq *coarse_seq;
    struct cbp_compressed_seq *cseq;
    struct cbp_seed_loc *seeds, *seedLoc;
    struct cbp_alignment alignment;
    char *kmer;
    int32_t seed_size, resind, mext, new_coarse_seq_id;
    int32_t last_match, current;
    int32_t id;
    bool has_end, changed;

    int32_t num_seeds = 0;

    cseq = cbp_compressed_seq_init(org_seq->id, org_seq->name);
    seed_size = coarse_db->seeds->seed_size;
    mext = compress_flags.match_extend;

    last_match = 0;
    current = 0;
    for (current = 0; current < org_seq->length - seed_size; current++) {
        kmer = org_seq->residues + current;
        seeds = cbp_seeds_lookup(coarse_db->seeds, kmer);
        if (seeds == NULL)
            continue;

        for (seedLoc = seeds; seedLoc != NULL; seedLoc = seedLoc->next) {
            num_seeds++;
            resind = seedLoc->residue_index;
            coarse_seq = cbp_coarse_get(coarse_db, seedLoc->coarse_seq_id);

            mlens = extend_match(
                mem,
                coarse_seq->seq->residues, resind, coarse_seq->seq->length,
                org_seq->residues, current, org_seq->length);

            has_end = mlens.olen + mext >= org_seq->length - current;
            if (mlens.olen < compress_flags.min_match_len && !has_end)
                continue;

            alignment = cbp_align_nw(
                mem,
                coarse_seq->seq->residues, resind, resind + mlens.rlen,
                org_seq->residues, current, current + mlens.olen);
            id = cbp_align_identity(
                alignment.ref, 0, alignment.length,
                alignment.org, 0, alignment.length);
            if (id < compress_flags.match_seq_id_threshold)
                continue;

            changed = false;
            if (has_end) {
                mlens.olen = org_seq->length - current;
                changed = true;
            }
            if (current - last_match <= mext) {
                mlens.olen += current - last_match;
                current = last_match;
                changed = true;
            }

            if (changed)
                alignment = cbp_align_nw(
                    mem,
                    coarse_seq->seq->residues, resind, resind + mlens.rlen,
                    org_seq->residues, current, current + mlens.olen);

            if (current - last_match > 0) {
                new_coarse_seq_id = add_without_match(
                    coarse_db, org_seq, last_match, current);
                cbp_compressed_seq_addlink(
                    cseq,
                    cbp_link_to_coarse_init_nodiff(
                        new_coarse_seq_id, 0, current - last_match));
            }

            cbp_compressed_seq_addlink(
                cseq,
                cbp_link_to_coarse_init(
                    coarse_seq->id, resind, resind + mlens.rlen, alignment));
            cbp_coarse_seq_addlink(
                coarse_seq,
                cbp_link_to_compressed_init(
                    org_seq->id, resind, resind + mlens.rlen));

            last_match = current + mlens.olen;
            current = last_match - 1;

            break;
        }

        cbp_seed_loc_free(seeds);
    }

    if (org_seq->length - last_match > 0) {
        new_coarse_seq_id = add_without_match(
            coarse_db, org_seq, last_match, org_seq->length);
        cbp_compressed_seq_addlink(
            cseq,
            cbp_link_to_coarse_init_nodiff(
                new_coarse_seq_id, 0, org_seq->length - last_match));
    }

    /* printf("Inspected %d seeds.\n", num_seeds); */

    return cseq;
}

struct extend_match
extend_match(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend,
             char *oseq, int32_t ostart, int32_t oend)
{
    struct cbp_alignment alignment;
    struct extend_match mlens;
    int32_t id;
    int32_t gwsize;
    int32_t rlen, olen;
    int32_t m;

    gwsize = compress_flags.gapped_window_size;
    rlen = rend - rstart;
    olen = oend - ostart;

    mlens.rlen = 0;
    mlens.olen = 0;
    while (true) {
        if (mlens.rlen == rlen || mlens.olen == olen)
            break;

        m = cbp_align_ungapped(
            compress_flags.ungapped_window_size,
            compress_flags.match_kmer_size,
            compress_flags.ext_seq_id_threshold,
            rseq, rstart + mlens.rlen, rend,
            oseq, ostart + mlens.olen, oend);

        mlens.rlen += m;
        mlens.olen += m;

        alignment = cbp_align_nw(
            mem,
            rseq, rstart + mlens.rlen, min(rend, rstart + mlens.rlen + gwsize),
            oseq, ostart + mlens.olen, min(oend, ostart + mlens.olen + gwsize));

        id = cbp_align_identity(
            alignment.ref, 0, alignment.length,
            alignment.org, 0, alignment.length);
        if (id < compress_flags.ext_seq_id_threshold)
            break;

        mlens.rlen += cbp_align_length_nogaps(alignment.ref);
        mlens.olen += cbp_align_length_nogaps(alignment.org);
    }

    return mlens;
}

static int32_t
add_without_match(struct cbp_coarse *coarse_db,
                  struct cbp_seq *org_seq, int32_t ostart, int32_t oend)
{
    struct cbp_coarse_seq *coarse_seq;

    coarse_seq = cbp_coarse_add(coarse_db, org_seq->residues, ostart, oend);
    cbp_coarse_seq_addlink(
        coarse_seq,
        cbp_link_to_compressed_init(org_seq->id, 0, oend - ostart));
    return coarse_seq->id;
}

static void *
cbp_compress_worker(void *data)
{
    struct worker_args *args;
    struct cbp_align_nw_memory *mem;
    struct cbp_seq *s;
    struct cbp_compressed_seq *cseq;

    args = (struct worker_args *) data;
    mem = cbp_align_nw_memory_init();
    while (NULL != (s = (struct cbp_seq *) ds_queue_pop(args->jobs))) {
        cseq = cbp_compress(args->db->coarse_db, s, mem);
        cbp_compressed_write(args->db->com_db, cseq);
        cbp_seq_free(s);
        cbp_compressed_seq_free(cseq);
    }
    cbp_align_nw_memory_free(mem);

    return NULL;
}

static int32_t
min(int32_t a, int32_t b)
{
    if (a < b)
        return a;
    return b;
}
