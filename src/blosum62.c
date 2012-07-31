#include <assert.h>
#include <stdbool.h>

#include "blosum62.h"

/* This is not good. It needs to be extremely fast, so anything like a hash
 * map is out of the question. The question is, how do we map a non-contiguous
 * alphabet to matrix indices?
 *
 * I'd also settle for this switch statement being auto-generated, but my
 * Ruby-fu is lacking. */
int32_t
blosum62_residue_to_index(char residue)
{
    switch (residue) {
    case 'A': return 0;
    case 'B': return 1;
    case 'C': return 2;
    case 'D': return 3;
    case 'E': return 4;
    case 'F': return 5;
    case 'G': return 6;
    case 'H': return 7;
    case 'I': return 8;
    /* case 'J': return ???; */
    case 'K': return 9;
    case 'L': return 10;
    case 'M': return 11;
    case 'N': return 12;
    /* case 'O': return ???; */
    case 'P': return 13;
    case 'Q': return 14;
    case 'R': return 15;
    case 'S': return 16;
    case 'T': return 17;
    /* case 'U': return ???; */
    case 'V': return 18;
    case 'W': return 19;
    case 'X': return 20;
    case 'Y': return 21;
    case 'Z': return 22;
    }

    assert(false);
    return 0;
}

