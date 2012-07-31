#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "opt.h"

#include "blosum62.h"
#include "database.h"
#include "flags.h"
#include "fasta.h"
#include "util.h"

#define FILENAME "../../cablastp/data/orf_trans_all.fasta"
/* #define FILENAME "/home/andrew/nr.fasta" */

struct job {
    struct fasta_seq_gen *fsg;
    char* id;
};

void * test(void *data);

int
main(int argc, char **argv)
{
    struct cbp_database *db;
    struct opt_config *conf;
    struct opt_args *args;
    int i;
    /* int cpus = 0; */
    /* struct fasta_seq_gen *fsg; */
    /* pthread_t *threads; */
    /* struct job *jobs; */

    conf = load_compress_args();
    args = opt_config_parse(conf, argc, argv);
    if (args->nargs < 2) {
        fprintf(stderr, 
            "Usage: %s [flags] database-dir fasta-file [ fasta-file ... ]\n",
            argv[0]);
        opt_config_print_usage(conf);
        exit(1);
    }

    db = cbp_database_init(args->args[0], compress_flags.seed_size, false);
    for (i = 1; i < args->nargs; i++)
        printf("Fasta file argument %d: %s\n", i, args->args[i]);

    /* fsg = fasta_generator_start(FILENAME, FASTA_EXCLUDE_NCBI_BLOSUM62, 10000); */
/*  */
    /* assert(threads = malloc(cpus * sizeof(*threads))); */
    /* assert(jobs = malloc(cpus * sizeof(*jobs))); */
    /* for (i = 0; i < num_cpus(); i++) { */
        /* char *id; */
/*  */
        /* assert(id = malloc(10 * sizeof(*id))); */
        /* sprintf(id, "%d", i); */
/*  */
        /* jobs[i].fsg = fsg; */
        /* jobs[i].id = id;  */
/*  */
        /* pthread_create(&(threads[i]), NULL, test, (void*) &(jobs[i])); */
    /* } */
    /* for (i = 0; i < cpus; i++) { */
        /* assert(0 == pthread_join(threads[i], NULL)); */
        /* free(jobs[i].id); */
    /* } */

    cbp_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);
    /* fasta_generator_free(fsg); */
    /* free(jobs); */
    /* free(threads); */

    return 0;
}

void *
test(void *data)
{
    struct fasta_seq *seq;
    struct job *j;
    int cnt = 0;
    unsigned int i;

    j = (struct job *) data;
    while (NULL != (seq = fasta_generator_next(j->fsg))) {
        cnt++;
        for (i = 0; i < strlen(seq->seq) * 2000; i++);
        fasta_free_seq(seq);
    }
    printf("Thread %s consumed %d sequences.\n", j->id, cnt);

    return NULL;
}
