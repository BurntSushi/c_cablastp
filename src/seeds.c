#include <stdio.h>
#include <stdlib.h>

#include "coarse.h"
#include "seeds.h"

const int8_t cbp_seeds_alpha_size[] = {
	0,  /* 'A' */
	1,  /* 'B' */
	2,  /* 'C' */
	3,  /* 'D' */
	4,  /* 'E' */
	5,  /* 'F' */
	6,  /* 'G' */
	7,  /* 'H' */
	8,  /* 'I' */
	-1, /* 'J' */
	9,  /* 'K' */
	10, /* 'L' */
	11, /* 'M' */
	12, /* 'N' */
	-1, /* 'O' */
	13, /* 'P' */
	14, /* 'Q' */
	15, /* 'R' */
	16, /* 'S' */
	17, /* 'T' */
	-1, /* 'U' */
	18, /* 'V' */
	19, /* 'W' */
	20, /* 'X' */
	21, /* 'Y' */
	22 /* 'Z' */
};

static int32_t residue_value(char residue);

static int32_t hash_kmer(struct cbp_seeds *seeds, char *kmer);

struct cbp_seeds *
cbp_seeds_init(int32_t seed_size)
{
    struct cbp_seeds *seeds;
    int32_t errno;
    int32_t i, p;

    seeds = malloc(sizeof(*seeds));
    assert(seeds);

    if (0 != (errno = pthread_rwlock_init(&seeds->lock, NULL))) {
        fprintf(stderr, "Could not create rwlock. Errno: %d\n", errno);
        exit(1);
    }

    seeds->seed_size = seed_size;
    seeds->powers_length = seed_size + 1;
    
    seeds->powers = malloc((seeds->powers_length) * sizeof(*seeds->powers));
    assert(seeds->powers);

    p = 1;
    for (i = 0; i < seeds->powers_length; i++) {
        seeds->powers[i] = p;
        p *= CABLASTP_SEEDS_ALPHA_SIZE;
    }

    seeds->locs_length = seeds->powers[seed_size];
    seeds->locs = malloc(seeds->locs_length * sizeof(*seeds->locs));
    assert(seeds->locs);
    for (i = 0; i < seeds->locs_length; i++)
        seeds->locs[i] = NULL;

    return seeds;
}

void
cbp_seeds_free(struct cbp_seeds *seeds) {
    int32_t errno;
    int32_t i;

    if (0 != (errno = pthread_rwlock_destroy(&seeds->lock))) {
        fprintf(stderr, "Could not destroy rwlock. Errno: %d\n", errno);
        exit(1);
    }
    for (i = 0; i < seeds->locs_length; i++)
        cbp_seed_loc_free(seeds->locs[i]); /* frees the whole list */
    free(seeds->locs);
    free(seeds->powers);
    free(seeds);
}

void
cbp_seeds_add(struct cbp_seeds *seeds, struct cbp_coarse_seq *seq)
{
    struct cbp_seed_loc *sl1, *sl2;
    char *kmer;
    int32_t hash, i;

    pthread_rwlock_wrlock(&seeds->lock);

    for (i = 0; i < seq->seq->length - seeds->seed_size; i++) {
        kmer = seq->seq->residues + i;
        sl1 = cbp_seed_loc_init(seq->id, i);

        hash = hash_kmer(seeds, kmer);
        if (seeds->locs[hash] == NULL)
            seeds->locs[hash] = sl1;
        else {
            for (sl2 = seeds->locs[hash]; sl2->next != NULL; sl2 = sl2->next);
            sl2->next = sl1;
        }
    }

    pthread_rwlock_unlock(&seeds->lock);
}

struct cbp_seed_loc *
cbp_seeds_lookup(struct cbp_seeds *seeds, char *kmer)
{
    struct cbp_seed_loc *sl, *copy_first, *copy;

    pthread_rwlock_rdlock(&seeds->lock);

    sl = seeds->locs[hash_kmer(seeds, kmer)];
    if (sl == NULL) {
        pthread_rwlock_unlock(&seeds->lock);
        return NULL;
    }

    copy_first = cbp_seed_loc_init(sl->coarse_seq_id, sl->residue_index);
    copy = copy_first;
    for (sl = sl->next; sl != NULL; sl = sl->next) {
        copy->next = cbp_seed_loc_init(sl->coarse_seq_id, sl->residue_index);
        copy = copy->next;
    }

    pthread_rwlock_unlock(&seeds->lock);

    return copy_first;
}

struct cbp_seed_loc *
cbp_seed_loc_init(uint32_t coarse_seq_id, uint16_t residue_index)
{
    struct cbp_seed_loc *seedLoc;

    seedLoc = malloc(sizeof(*seedLoc));
    assert(seedLoc);

    seedLoc->coarse_seq_id = coarse_seq_id;
    seedLoc->residue_index = residue_index;
    seedLoc->next = NULL;

    return seedLoc;
}

void
cbp_seed_loc_free(struct cbp_seed_loc *seedLoc)
{
    struct cbp_seed_loc *seed1, *seed2;

    for (seed1 = seedLoc; seed1 != NULL; ) {
        seed2 = seed1->next;
        free(seed1);
        seed1 = seed2;
    }
}

static int32_t residue_value(char residue)
{
    int32_t i, val;

    i = residue - 'A';
    if (i < 0 || i >= 26) {
        fprintf(stderr, "Invalid amino acid residue: %c\n", residue);
        exit(1);
    }
    if (-1 == (val = cbp_seeds_alpha_size[i])) {
        fprintf(stderr, "Invalid amino acid residue: %c\n", residue);
        exit(1);
    }
    return val;
}

static int32_t hash_kmer(struct cbp_seeds *seeds, char *kmer)
{
    int32_t i, key;

    key = 0;
    for (i = 0; i < seeds->seed_size; i++)
        key += residue_value(kmer[i]) * seeds->powers[i];
    return key;
}
