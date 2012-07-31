#ifndef __CABLASTP_SEEDS_H__
#define __CABLASTP_SEEDS_H__

/* Apparently this is required to make pthread_rwlock* stuff available. */
#define __USE_UNIX98

#include <pthread.h>
#include <stdint.h>

#define CABLASTP_SEEDS_ALPHA_SIZE 23

const int8_t cbp_seeds_alpha_size[26];

struct cbp_seed_loc {
    uint32_t coarse_seq_id;
    uint16_t residue_index;
    struct cbp_seed_loc *next;
};

struct cbp_seed_loc *
cbp_seed_loc_init(uint32_t coarse_seq_id, uint16_t residue_index);

void
cbp_seed_loc_free(struct cbp_seed_loc *seedLoc);

struct cbp_seeds {
    int32_t seed_size;
    struct cbp_seed_loc **locs;
    int32_t locs_length;
    int32_t *powers;
    int32_t powers_length;
    pthread_rwlock_t lock;
};

struct cbp_seeds *
cbp_seeds_init(int32_t seed_size);

void
cbp_seeds_free(struct cbp_seeds *seeds);

struct cbp_coarse_seq;

void
cbp_seeds_add(struct cbp_seeds *seeds, struct cbp_coarse_seq *seq);

/* Produces a copy of the list of seeds for 'kmer', and therefore the
 * result needs to be freed with `cbp_seed_loc_free` when finished. */
struct cbp_seed_loc *
cbp_seeds_lookup(struct cbp_seeds *seeds, char *kmer);

#endif
