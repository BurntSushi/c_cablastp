#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "ds.h"

#include "compressed.h"

static struct cbp_compressed_seq *
seq_at(struct cbp_compressed *com_db, int32_t i);

struct cbp_compressed *
cbp_compressed_init(FILE *file_compressed, FILE *file_index)
{
    struct cbp_compressed *com_db;

    com_db = malloc(sizeof(*com_db));
    assert(com_db);

    com_db->file_compressed = file_compressed;
    com_db->file_index = file_index;
    com_db->seqs = ds_vector_create_capacity(100);

    return com_db;
}

void
cbp_compressed_free(struct cbp_compressed *com_db)
{
    int i;

    fclose(com_db->file_compressed);
    fclose(com_db->file_index);

    for (i = 0; i < com_db->seqs->size; i++)
        cbp_compressed_seq_free(seq_at(com_db, i));

    ds_vector_free_no_data(com_db->seqs);
    free(com_db);
}

int32_t
cbp_compressed_size(struct cbp_compressed *com_db)
{
    return com_db->seqs->size;
}

void
cbp_compressed_add(struct cbp_compressed *com_db,
                   struct cbp_compressed_seq *seq)
{
    ds_vector_append(com_db->seqs, (void*) seq);
}

void
cbp_compressed_save_plain(struct cbp_compressed *com_db)
{
    int i;
    struct cbp_compressed_seq *seq;
    struct cbp_link_to_coarse *link;

    for (i = 0; i < com_db->seqs->size; i++) {
        seq = seq_at(com_db, i);
        fprintf(com_db->file_compressed, "> %d; %s\n", seq->id, seq->name);
        for (link = seq->links; link != NULL; link = link->next)
            fprintf(com_db->file_compressed,
                "reference sequence id: %d, reference range: (%d, %d)\n%s",
                link->coarse_seq_id, link->coarse_start, link->coarse_end,
                link->diff);
    }
}

struct cbp_compressed_seq *
cbp_compressed_seq_init(int32_t id, char *name)
{
    struct cbp_compressed_seq *seq;

    seq = malloc(sizeof(*seq));
    assert(seq);

    seq->id = id;
    seq->links = NULL;

    seq->name = malloc((1 + strlen(name)) * sizeof(*seq->name));
    assert(seq->name);
    strcpy(seq->name, name);

    return seq;
}

void
cbp_compressed_seq_free(struct cbp_compressed_seq *seq)
{
    struct cbp_link_to_coarse *link1, *link2;

    for (link1 = seq->links; link1 != NULL; ) {
        link2 = link1->next;
        cbp_link_to_coarse_free(link1);
        link1 = link2;
    }
    free(seq->name);
    free(seq);
}

void
cbp_compressed_seq_addlink(struct cbp_compressed_seq *seq,
                           struct cbp_link_to_coarse *newlink)
{
    struct cbp_link_to_coarse *link;

    assert(newlink->next == NULL);

    if (seq->links == NULL) {
        seq->links = newlink;
        return;
    }

    for (link = seq->links; link->next != NULL; link = link->next);
    link->next = newlink;
}


struct cbp_link_to_coarse *
cbp_link_to_coarse_init(int32_t coarse_seq_id,
                        int16_t coarse_start, int16_t coarse_end,
                        char *align_ref, char *align_tgt) {
    struct cbp_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->diff = NULL;
    link->coarse_seq_id = coarse_seq_id;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;

    return link;
}

struct cbp_link_to_coarse *
cbp_link_to_coarse_init_nodiff(int32_t coarse_seq_id,
                               int16_t coarse_start, int16_t coarse_end) {
    struct cbp_link_to_coarse *link;

    link = malloc(sizeof(*link));
    assert(link);

    link->diff = NULL;
    link->coarse_seq_id = coarse_seq_id;
    link->coarse_start = coarse_start;
    link->coarse_end = coarse_end;
    link->next = NULL;

    return link;
}

void
cbp_link_to_coarse_free(struct cbp_link_to_coarse *link)
{
    free(link->diff);
    free(link);
}

static struct cbp_compressed_seq *
seq_at(struct cbp_compressed *com_db, int32_t i)
{
    return (struct cbp_compressed_seq *) ds_vector_get(com_db->seqs, i);
}
