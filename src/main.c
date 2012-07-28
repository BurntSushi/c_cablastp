#include <stdio.h>

#include "blosum62.h"
#include "fasta.h"

int
/* main(int argc, char **argv) */
main(void)
{
    int i;
    struct fasta_file *ff;

    ff = fasta_read_all("../cablastp/data/small.fasta",
        FASTA_EXCLUDE_NCBI_BLOSUM62);
    for (i = 0; i < ff->length; i++) {
        printf("%s\n%s\n", ff->seqs[0]->name, ff->seqs[0]->seq);
    }

    fasta_free_all(ff);
    return 0;
}
