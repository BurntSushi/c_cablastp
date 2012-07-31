#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "database.h"

static FILE * open_db_file(char *path);

static char * path_join(char *a, char *b);

static char * basename(char *path);

struct cbp_database *
cbp_database_init(char *dir, int32_t seed_size, bool add)
{
    struct cbp_database *db;
    struct stat buf;
    FILE *ffasta, *fseeds, *flinks, *fcompressed, *findex;
    char *pfasta, *pseeds, *plinks, *pcompressed, *pindex;

    /* If we're not adding to a database, make sure `dir` does not exist. */
    if (!add && 0 == stat(dir, &buf)) {
        fprintf(stderr,
            "The directory '%s' already exists. A new compressed "
            "database cannot be created in the same directory as an "
            "existing database. If you want to append to an existing "
            "database, use the '--append' flag.\n", dir);
        exit(1);
    }
    /* Otherwise, check to make sure it *does* exist. */
    if (add && 0 != stat(dir, &buf)) {
        fprintf(stderr, "Could not open '%s' database for appending.", dir);
        exit(1);
    }

    if (0 != mkdir(dir, 0777)) {
        fprintf(stderr, "cbp_database_init: 'mkdir %s' failed because: %s\n",
            dir, strerror(errno));
        exit(1);
    }

    db = malloc(sizeof(*db));
    assert(db);
    db->name = basename(dir);

    pfasta = path_join(dir, CABLASTP_COARSE_FASTA);
    pseeds = path_join(dir, CABLASTP_COARSE_LINKS);
    plinks = path_join(dir, CABLASTP_COARSE_SEEDS);
    pcompressed = path_join(dir, CABLASTP_COMPRESSED);
    pindex = path_join(dir, CABLASTP_INDEX);

    ffasta = open_db_file(pfasta);
    fseeds = open_db_file(pseeds);
    flinks = open_db_file(plinks);
    fcompressed = open_db_file(pcompressed);
    findex = open_db_file(pindex);

    db->coarse_db = cbp_coarse_init(seed_size, ffasta, fseeds, flinks);
    db->com_db = cbp_compressed_init(fcompressed, findex);

    free(pfasta);
    free(pseeds);
    free(plinks);
    free(pcompressed);
    free(pindex);

    return db;
}

void
cbp_database_free(struct cbp_database *db)
{
    /* All files opened in cbp_database_init are close in subsequent frees. */
    cbp_coarse_free(db->coarse_db);
    cbp_compressed_free(db->com_db);
    free(db->name);
    free(db);
}

static FILE *
open_db_file(char *path)
{
    struct stat buf;
    FILE *fp;

    if (0 != stat(path, &buf))
        if (-1 == creat(path, 0666)) {
            fprintf(stderr, "open_db_file: 'creat %s' failed: %s\n",
                path, strerror(errno));
            exit(1);
        }
    if (NULL == (fp = fopen(path, "r+"))) {
        fprintf(stderr, "open_db_file: 'fopen %s' failed: %s\n",
            path, strerror(errno));
        exit(1);
    }

    return fp;
}

static char *
path_join(char *a, char *b)
{
    char *joined;

    joined = malloc((1 + strlen(a) + 1 + strlen(b)) * sizeof(*joined));
    assert(joined);

    sprintf(joined, "%s/%s", a, b);
    return joined;
}

static char *
basename(char *path)
{
    char *base;
    int i;
    int len;

    len = strlen(path);
    for (i = len; i >= 0 && path[i] != '/'; i--);
    if (i > 0)
        i++;

    base = malloc((1 + len - i) * sizeof(*base));
    assert(base);

    strncpy(base, path + i, len - i);

    return base;
}
