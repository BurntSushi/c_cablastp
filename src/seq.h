#ifndef __CABLASTP_SEQ_H__
#define __CABLASTP_SEQ_H__

#include <stdint.h>

struct cbp_seq {
    int32_t id;
    char *name;
    char *residues;
    int32_t length;
};

struct cbp_seq *
cbp_seq_init(int32_t id, char *name, char *residues);

struct cbp_seq *
cbp_seq_init_range(int32_t id, char *name, char *residues,
                   int32_t start, int32_t end);

void
cbp_seq_free(struct cbp_seq *seq);

#endif
