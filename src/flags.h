#ifndef __CABLASTP_OPT_H__
#define __CABLASTP_OPT_H__

#include <stdint.h>

#include "opt.h"

struct compress_flags {
    int32_t gapped_window_size;
    int32_t ungapped_window_size;
    int32_t match_extend;
    int32_t match_kmer_size;
    int32_t match_seq_id_threshold;
    int32_t min_match_len;
    int32_t procs;
    int32_t seed_size;
    int32_t ext_seq_id_threshold;
} compress_flags;

struct opt_config *
load_compress_args();

#endif
