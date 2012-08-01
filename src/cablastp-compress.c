#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

#include "opt.h"

#include "blosum62.h"
#include "compression.h"
#include "database.h"
#include "flags.h"
#include "fasta.h"
#include "seq.h"
#include "util.h"

int
main(int argc, char **argv)
{
    struct cbp_database *db;
    struct cbp_compress_workers *workers;
    struct opt_config *conf;
    struct opt_args *args;
    struct fasta_seq_gen *fsg;
    struct fasta_seq *seq;
    struct cbp_seq *org_seq;
    int i, org_seq_id;
    struct timeval start, current;
    long double elapsed;

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
    workers = cbp_compress_start_workers(db, compress_flags.procs);

    org_seq_id = 0;
    gettimeofday(&start, NULL);
    for (i = 1; i < args->nargs; i++) {
        fsg = fasta_generator_start(
            args->args[i], FASTA_EXCLUDE_NCBI_BLOSUM62, 100);

        while (NULL != (seq = fasta_generator_next(fsg))) {
            org_seq = cbp_seq_init(org_seq_id, seq->name, seq->seq);
            cbp_compress_send_job(workers, org_seq);

            fasta_free_seq(seq);

            org_seq_id++;
            if (org_seq_id % 1000 == 0) {
                gettimeofday(&current, NULL);
                elapsed = (long double)(current.tv_sec - start.tv_sec);
                printf("%d sequences compressed (%0.4Lf seqs/sec)\n",
                    org_seq_id, ((long double) org_seq_id) / elapsed);
            }
        }

        fasta_generator_free(fsg);
    }

    cbp_compress_join_workers(workers);
    cbp_coarse_save_plain(db->coarse_db);
    cbp_compressed_save_plain(db->com_db);

    cbp_database_free(db);
    opt_config_free(conf);
    opt_args_free(args);
    cbp_compress_free_workers(workers);

    return 0;
}
