#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "blosum62.h"
#include "fasta.h"

#define FILENAME "../cablastp/data/orf_trans_all.fasta"
/* #define FILENAME "/home/andrew/nr.fasta" */

struct job {
    struct fasta_seq_gen *fsg;
    char* id;
};

void * test(void *data);

int
/* main(int argc, char **argv) */
main(void)
{
    struct fasta_seq_gen *fsg;
    pthread_t *threads;
    struct job *jobs;
    int i;
    int cpus = 0;

    cpus = (int) sysconf(_SC_NPROCESSORS_ONLN);
    fsg = fasta_generator_start(FILENAME, FASTA_EXCLUDE_NCBI_BLOSUM62, 1000);

    assert(threads = malloc(cpus * sizeof(*threads)));
    assert(jobs = malloc(cpus * sizeof(*jobs)));
    for (i = 0; i < cpus; i++) {
        char *id;

        assert(id = malloc(10 * sizeof(*id)));
        sprintf(id, "%d", i);

        jobs[i].fsg = fsg;
        jobs[i].id = id; 

        pthread_create(&(threads[i]), NULL, test, (void*) &(jobs[i]));
    }
    for (i = 0; i < cpus; i++) {
        assert(0 == pthread_join(threads[i], NULL));
        free(jobs[i].id);
    }

    fasta_generator_free(fsg);
    free(jobs);
    free(threads);

    return 0;
}

void *
test(void *data)
{
    struct fasta_seq *seq;
    struct job *j;
    int cnt = 0;
    /* unsigned int i; */

    j = (struct job *) data;
    while (NULL != (seq = fasta_generator_next(j->fsg))) {
        cnt++;
        /* for (i = 0; i < strlen(seq->seq) * 2000; i++); */
        fasta_free_seq(seq);
    }
    printf("Thread %s consumed %d sequences.\n", j->id, cnt);

    return NULL;
}
