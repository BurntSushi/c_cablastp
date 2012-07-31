#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "coarse.h"

struct cbp_coarse *
cbp_coarse_init(int32_t seed_size,
                FILE *file_fasta, FILE *file_seeds, FILE *file_links)
{
    struct cbp_coarse *coarse_db;
    int32_t errno;

    coarse_db = malloc(sizeof(*coarse_db));
    assert(coarse_db);

    coarse_db->seqs = ds_vector_create_capacity(1000000);
    /* coarse_db->seeds = cbp_seeds_init(seed_size); */
    coarse_db->file_fasta = file_fasta;
    coarse_db->file_seeds = file_seeds;
    coarse_db->file_links = file_links;

    if (0 != (errno = pthread_rwlock_init(&coarse_db->lock_seq, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return coarse_db;
}

void
cbp_coarse_free(struct cbp_coarse *coarse_db)
{
    int32_t errno;
    int32_t i;

    fclose(coarse_db->file_fasta);
    fclose(coarse_db->file_seeds);
    fclose(coarse_db->file_links);

    if (0 != (errno = pthread_rwlock_destroy(&coarse_db->lock_seq))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (i = 0; i < coarse_db->seqs->size; i++)
        cbp_coarse_seq_free(
            (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i));

    ds_vector_free_no_data(coarse_db->seqs);
    free(coarse_db);
}

struct cbp_coarse_seq *
cbp_coarse_add(struct cbp_coarse *coarse_db,
               char *residues, int32_t start, int32_t end)
{
    struct cbp_coarse_seq *seq;
    int32_t id;

    pthread_rwlock_wrlock(&coarse_db->lock_seq);
    id = coarse_db->seqs->size;
    seq = cbp_coarse_seq_init(id, residues, start, end);
    ds_vector_append(coarse_db->seqs, (void*) seq);
    pthread_rwlock_unlock(&coarse_db->lock_seq);

    /* add to seeds */

    return seq;
}

struct cbp_coarse_seq *
cbp_coarse_get(struct cbp_coarse *coarse_db, int32_t i)
{
    struct cbp_coarse_seq *seq;

    pthread_rwlock_rdlock(&coarse_db->lock_seq);
    seq = (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
    pthread_rwlock_unlock(&coarse_db->lock_seq);

    return seq;
}

void
cbp_coarse_save_plain(struct cbp_coarse *coarse_db)
{
    struct cbp_coarse_seq *seq;
    struct cbp_link_to_compressed *link;
    int32_t i;

    for (i = 0; i < coarse_db->seqs->size; i++) {
        seq = (struct cbp_coarse_seq *) ds_vector_get(coarse_db->seqs, i);
        fprintf(coarse_db->file_fasta, "> %d\n%s\n", i, seq->seq->residues);

        fprintf(coarse_db->file_links, "> %d\n", i);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(coarse_db->file_links,
                "original sequence id: %d, reference range: (%d, %d)\n",
                link->org_seq_id, link->coarse_start, link->coarse_end);
    }
}

struct cbp_coarse_seq *
cbp_coarse_seq_init(int32_t id, char *residues, int32_t start, int32_t end)
{
    struct cbp_coarse_seq *seq;
    int32_t errno;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->seq = cbp_seq_init_range(id, "", residues, start, end);
    seq->links = NULL;

    if (0 != (errno = pthread_rwlock_init(&seq->lock_links, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    return seq;
}

void
cbp_coarse_seq_free(struct cbp_coarse_seq *seq)
{
    struct cbp_link_to_compressed *link1, *link2;
    int32_t errno;

    if (0 != (errno = pthread_rwlock_destroy(&seq->lock_links))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (link1 = seq->links; link1 != NULL; ) {
        link2 = link1->next;
        cbp_link_to_compressed_free(link1);
        link1 = link2;
    }

    cbp_seq_free(seq->seq);
    free(seq);
}

void
cbp_coarse_seq_addlink(struct cbp_coarse_seq *seq,
                       struct cbp_link_to_compressed *newlink)
{
    struct cbp_link_to_compressed *link;

    assert(newlink->next == NULL);

    pthread_rwlock_wrlock(&seq->lock_links);
    if (seq->links == NULL) {
        seq->links = newlink;
        pthread_rwlock_unlock(&seq->lock_links);
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
    pthread_rwlock_unlock(&seq->lock_links);
}

struct cbp_link_to_compressed *
cbp_link_to_compressed_init(int32_t org_seq_id,
                            int16_t coarse_start, int16_t coarse_end)
{
    struct cbp_link_to_compressed *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->org_seq_id = org_seq_id;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;

    return link;
}

void
cbp_link_to_compressed_free(struct cbp_link_to_compressed *link)
{
    free(link);
}
