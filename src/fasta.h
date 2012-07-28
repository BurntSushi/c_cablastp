#ifndef __CABLASTP_FASTA_H__
#define __CABLASTP_FASTA_H__

#include <stdio.h>
#include <stdint.h>

#define FASTA_INITIAL_SEQUENCE_LENGTH 1000
#define FASTA_MAX_LINE 1024

struct fasta_file {
    struct fasta_seq **seqs;
    int32_t length;
};

struct fasta_seq {
    char *name;
    char *seq;
};

struct fasta_file *
fasta_read_all(const char * file_name);

struct fasta_seq *
fasta_read_next(FILE *f);

void
fasta_free_all_seqs(struct fasta_file *ff);

void
fasta_free_seq(struct fasta_seq *seq);

#endif
