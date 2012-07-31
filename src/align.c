#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"
#include "blosum62.h"

int32_t
cbp_align_ungapped(int32_t window_size, int32_t kmer_size, int32_t id_threshold,
                   char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend)
{
    int32_t length, scanned, successive;
    int32_t rlen, olen;
    int32_t id;
    int32_t i, rs, os;
    bool try_next_window;

    length = 0;
    scanned = 0;
    successive = 0;
    rlen = rend - rstart;
    olen = oend - ostart;
    try_next_window = true;

    while (try_next_window) {
        try_next_window = false;
        for (i = 0; i < window_size; i++) {
            if (scanned >= rlen || scanned >= olen)
                break;

            if (rseq[rstart + scanned] == oseq[ostart + scanned])
                successive++;
            else
                successive = 0;

            scanned++;
            if (successive == kmer_size) {
                if ((scanned - kmer_size) - length > 0) {
                    rs = rstart + length;
                    os = ostart + length;
                    id = cbp_align_identity(
                        rseq, rs, rstart + scanned - kmer_size,
                        oseq, os, ostart + scanned - kmer_size);
                    if (id < id_threshold) {
                        successive--;
                        continue;
                    }
                }

                length = scanned;
                successive = 0;
                try_next_window = true;
                break;
            }
        }
    }
    return length;
}

int32_t
cbp_align_identity(char *rseq, int32_t rstart, int32_t rend,
                   char *oseq, int32_t ostart, int32_t oend)
{
    int32_t rlen, olen, same, i;

    rlen = rend - rstart;
    olen = oend - ostart;
    assert(rlen == olen);

    if (rlen == 0 && olen == 0)
        return 0;

    same = 0;
    for (i = 0; i < rlen; i++)
        if (rseq[rstart + i] == oseq[ostart + i])
            same++;

    return (same * 100) / rlen;
}

struct cbp_align_nw_memory *
cbp_align_nw_memory_init()
{
    struct cbp_align_nw_memory *mem;
    int seq_size = CABLASTP_ALIGN_SEQ_SIZE;

    mem = malloc(sizeof(*mem));
    assert(mem);

    mem->table = malloc(seq_size * seq_size * sizeof(*mem->table));
    assert(mem->table);

    mem->ref = malloc(seq_size * sizeof(*mem->ref));
    assert(mem->ref);

    mem->org = malloc(seq_size * sizeof(*mem->org));
    assert(mem->org);

    return mem;
}

void
cbp_align_nw_memory_free(struct cbp_align_nw_memory *mem)
{
    free(mem->table);
    free(mem->ref);
    free(mem->org);
    free(mem);
}

struct cbp_alignment
cbp_align_nw(struct cbp_align_nw_memory *mem,
             char *rseq, int32_t rstart, int32_t rend,
             char *oseq, int32_t ostart, int32_t oend)
{
    struct cbp_alignment alignment;
    int32_t *table;
    int32_t gapi;
    int32_t i, i2, i3, j, r, c, off, tablen;
    int32_t sdiag, sup, sleft, rval, oval;
    int32_t constraint;
    bool constrained;
    bool allocated = false;
    char tmp;

    gapi = BLOSUM62_SIZE - 1;
    r = rend - rstart + 1; /* include gap penalty */
    c = oend - ostart + 1; /* include gap penalty */
    tablen = r * c;
    off = 0;

    if (tablen > CABLASTP_ALIGN_SEQ_SIZE * CABLASTP_ALIGN_SEQ_SIZE) {
        table = malloc(tablen * sizeof(*table));
        assert(table);
        allocated = true;
    } else {
        table = mem->table;
        for (i = 0; i < tablen; i++)
            table[i] = 0;
    }

    constrained = true;
    constraint = r / 4;
    if (r <= 11 || c <= 11)
        constrained = false;

    for (i = 1; i < r; i++) {
        i2 = (i - 1) * c;
        i3 = i * c;
        for (j = 1; j < c; j++) {
            if (constrained && ((i-j) > constraint || (j-i) > constraint))
                continue;
            rval = blosum62_residue_to_index(rseq[rstart + i - 1]);
            oval = blosum62_residue_to_index(oseq[ostart + j - 1]);

            off = i2 + (j - 1);
            sdiag = table[off] + BLOSUM62_MATRIX[rval][oval];
            sup = table[off+1] + BLOSUM62_MATRIX[rval][gapi];
            sleft = table[off+c] + BLOSUM62_MATRIX[gapi][oval];
            
            if (sdiag > sup && sdiag > sleft)
                table[i3+j] = sdiag;
            else if (sup > sleft)
                table[i3+j] = sup;
            else
                table[i3+j] = sleft;
        }
    }

    alignment.ref = mem->ref;
    alignment.org = mem->org;
    alignment.length = 0;

    i = r - 1;
    j = c - 1;
    while (i > 0 && j > 0) {
        rval = blosum62_residue_to_index(rseq[rstart + i - 1]);
        oval = blosum62_residue_to_index(oseq[ostart + j - 1]);

        sdiag = table[(i-1) * c + (j-1)] + BLOSUM62_MATRIX[rval][oval];
        sup = table[(i-1) * c + j] + BLOSUM62_MATRIX[gapi][oval];
        sleft = table[i * c + (j-1)] + BLOSUM62_MATRIX[rval][gapi];

        if (sdiag > sup && sdiag > sleft) {
            i--;
            j--;
            alignment.ref[alignment.length] = rseq[rstart + i];
            alignment.org[alignment.length] = oseq[ostart + j];
            alignment.length++;
        } else if (sup > sleft) {
            i--;
            alignment.ref[alignment.length] = rseq[rstart + i];
            alignment.org[alignment.length] = '-';
            alignment.length++;
        } else {
            j--;
            alignment.ref[alignment.length] = '-';
            alignment.org[alignment.length] = oseq[ostart + j];
            alignment.length++;
        }
    }
    for ( ; i > 0; i--) {
        alignment.ref[alignment.length] = rseq[rstart + i - 1];
        alignment.org[alignment.length] = '-';
        alignment.length++;
    }
    for ( ; j > 0; j--) {
        alignment.ref[alignment.length] = '-';
        alignment.org[alignment.length] = oseq[ostart + j - 1];
        alignment.length++;
    }

    for (i = 0, j = alignment.length - 1; i < j; i++, j--) {
        tmp = alignment.ref[i];
        alignment.ref[i] = alignment.ref[j];
        alignment.ref[j] = tmp;

        tmp = alignment.org[i];
        alignment.org[i] = alignment.org[j];
        alignment.org[j] = tmp;
    }

    if (allocated)
        free(table);

    alignment.ref[alignment.length] = '\0';
    alignment.org[alignment.length] = '\0';

    return alignment;
}

int32_t
cbp_align_length_nogaps(char *residues)
{
    int i, len, rlen;

    len = 0;
    rlen = strlen(residues);
    for (i = 0; i < rlen; i++)
        if (residues[i] != '-')
            len++;

    return len;
}
