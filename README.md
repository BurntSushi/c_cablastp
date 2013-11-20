This is an early implementation of CaBLASTP in ANSI C. Its purpose was to 
roughly compare the performance of the Go implementation with an equivalent C 
implementation. Namely, this implementation roughly corresponds to commit
[219c79](https://github.com/BurntSushi/cablastp/commit/219c792c0ab2d5d46d791dbd237f4d09ea3fccf3)
in the official Go implementation.

The compressed databases created by this code are not in any way compatible 
with the compressed databases produced by the official implementation. Namely, 
the data itself should be the same (or similar), but its representation on disk 
is completely different. There may be other, subtler differences as well.

**This code should not be used.** It is provided for those interested in an
alternative implementation.

There are no benchmarks.

The official Go implementation can be found here:
https://github.com/BurntSushi/cablastp


Installation
============
c_cablastp depends on three libraries: `opt`, `ds` and `pthread`. `pthread` 
should be installed via your system's package manager. `opt` and `ds` can be 
found in my [clibs respository](https://github.com/BurntSushi/clibs).

Briefly, the following commands should get c_cablastp into a working state:

```bash
mkdir c_cablastp
git clone git://github.com/BurntSushi/clibs
git clone git://github.com/BurntSushi/c_cablastp
cd clibs
make
export C_INCLUDE_PATH=$(pwd)/include
export LIBRARY_PATH=$(pwd)/lib
cd ../c_cablastp
make
./cablastp-compress --help
```

The usage of `cablastp-compress` is:

```bash
cablastp-compress [flags] database-directory fasta-file [fasta-file ...]
```

