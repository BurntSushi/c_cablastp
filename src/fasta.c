#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"

static bool
is_new_sequence_start(FILE *f);

static void
chomp(char *s1);

struct fasta_file *
fasta_read_all(const char * file_name)
{
    FILE *f;
    struct fasta_file *ff;
    struct fasta_seq *new_seq;
    int allocated = 0;

    printf("Reading all sequences from: %s\n", file_name);

    f = fopen(file_name, "r");
    if (f == NULL) {
        perror("fasta_read_all");
        exit(1);
    }

    ff = malloc(sizeof(*ff));
    assert(ff);

    allocated = FASTA_INITIAL_SEQUENCE_LENGTH;
    ff->seqs = malloc(allocated * sizeof(*ff->seqs));
    assert(ff->seqs);

    while (NULL != (new_seq = fasta_read_next(f))) {
        if (ff->length == allocated) {
            allocated *= 1.5;
            ff->seqs = realloc(ff->seqs, allocated * sizeof(*ff->seqs));
            assert(ff->seqs);
        }
        ff->seqs[ff->length++] = new_seq;
    }
    
    return ff;
}

struct fasta_seq *
fasta_read_next(FILE *f)
{
    struct fasta_seq *fs;
    char buf[FASTA_MAX_LINE];
    char *fgets_status;

    fs = malloc(sizeof(*fs));
    assert(fs);

    /* check to make sure the next line starts a new sequence record */
    if (!is_new_sequence_start(f))
        return NULL;

    /* read in the sequence id */
    fgets_status = fgets(buf, FASTA_MAX_LINE, f);
    if (fgets_status == NULL)
        return NULL;

    chomp(buf);

    fs->name = malloc((strlen(buf) + 1) * sizeof(*fs->name));
    assert(fs->name);
    strcpy(fs->name, buf);

    /* Now read all of the sequence data for this record */
    fs->seq = malloc(sizeof(*fs->seq));
    assert(fs->seq);
    fs->seq[0] = '\0';
    while (!is_new_sequence_start(f)) {
        memset(buf, 0, FASTA_MAX_LINE);

        fgets_status = fgets(buf, FASTA_MAX_LINE, f);
        if (fgets_status == NULL) {
            return fs;
        }

        chomp(buf);
        fs->seq = realloc(
            fs->seq, sizeof(*fs->seq) * (1 + strlen(buf) + strlen(fs->seq)));
        assert(fs->seq);
        strcat(fs->seq, buf);
    }

    return fs;
}

void
fasta_free_all_seqs(struct fasta_file *ff)
{
    int i;

    for (i = 0; i < ff->length; i++)
        fasta_free_seq(ff->seqs[i]);
    free(ff->seqs);
    free(ff);
}

void
fasta_free_seq(struct fasta_seq *seq)
{
    free(seq->name);
    free(seq->seq);
    free(seq);
}

static bool
is_new_sequence_start(FILE *f)
{
    char next;
    bool is_new_seq;

    next = fgetc(f);
    is_new_seq = next == '>';
    if (EOF == ungetc(next, f)) {
        return false;
        perror("is_new_sequence_start");
        exit(1);
    }

    return is_new_seq;
}

static void
chomp(char *s1)
{
    int len;

    len = strlen(s1);
    if (s1[len - 1] == '\n')
        s1[len - 1] = '\0';
    if (s1[len - 1] == '\r')
        s1[len - 1] = '\0';
}
