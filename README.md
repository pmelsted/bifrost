# BFGraph

Highly Parallel and Memory Efficient Compacted de Bruijn Graph Construction

This repository contains the source code for a new parallel and memory efficient algorithm enabling the direct construction of the compacted de Bruijn graph without producing the intermediate uncompacted de Bruijn graph. Despite making extensive use of a probabilistic data structure, the Bloom filter, our algorithm guarantees that the produced compacted de Bruijn graph is deterministic. Furthermore, the algorithm features de Bruijn graph simplification steps used by assemblers such as tip clipping and isolated unitig removal. In addition, as disk-based software performance is significantly affected by the discrepancy of speed among disk storage technologies, our method uses only main memory storage.

## Dependencies

In order to compile and use the BFT, you need a machine running a 64 bits Linux or Mac OS operating system. BFGraph successfully compiles and runs on Ubuntu 17.04 and MacOS.

In order to install BFGraph, you will need Cmake (https://cmake.org/), Jemalloc (http://www.canonware.com/jemalloc) and zlib (https://zlib.net/). All can be downloaded and installed by following the instructions on their respective websites. It is however most likely that at least few of them are available via a package manager for your operating system.

If you operating system is Ubuntu/Debian:
```
sudo apt-get install cmake libjemalloc1 libjemalloc-dev zlib1g
```

If you operating system is MacOS, Jemalloc can be easily downloaded and installed via Homebrew:
```
brew install cmake jemalloc zlib
```

## Compilation and Installation

Compilation requires Cmake (version 2.8.12 minimum) and a compiler such as GCC or Clang. You can verify the presence of GCC or Clang on your system with:
```
gcc -v
g++ -v
```

If not present (unlikely), they can be installed for Ubuntu/Debian with:
```
sudo apt-get install build-essential
```

or for MacOS:
```
brew install --with-toolchain llvm
```

Then, installing should be as easy as:
```
cd <BFGraph_directory>
mkdir build
cd build
cmake ..
make
```

The default maximum *k*-mer size supported is 31. To work with larger *k*, you must replace *MAX_KMER_SIZE* in *CMakeLists.txt* with an larger (appropriate) number, such as:
```
set( MAX_KMER_SIZE "64")
```

In this case, the maximum *k* allowed is 63. Keep in mind that increasing MAX_KMER_SIZE increases BFGraph memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

## Binary usage:

Type the following command to run BFGraph:
```
./BFGraph
```

It should display the list of available commands:
```
A memory efficient de Bruijn graph assembler.

Usage: BFGraph <cmd> [options] ...

Where <cmd> can be one of:
    filter       Filters errors from reads
    contigs      Builds an initial contig graph
    cite         Prints information for citing the paper
    version      Displays version number
```

BFGraph works (for now) in two steps: first, errors are removed from the reads (*filter* command), then, the compacted de Bruijn graph is built from the filtered reads (*contigs* command).

### Filtering

```
./BFGraph filter
```

```
Filter errors from FASTA/FASTQ files

Usage: BFGraph filter [arguments] ... FASTQ_files

Required arguments:
  -t, --threads=INT           Number of threads (default is 1)
  -k, --kmer-size=INT         Size of k-mers (default is 31)
  -g, --min-size=INT          Size of minimizers (default is 23)
      --ref                   Reference mode, no filtering
  -n, --num-kmers=LONG        Estimated number (upper bound) of different k-mers in the reads
  -b, --bloom-bits=INT        Number of Bloom filter bits per k-mer occurring at least once in the reads (default is 14)
  -N, --num-kmer2=LONG        Estimated number (upper bound) of different k-mers occurring at least twice in the reads
  -B, --bloom-bits2=INT       Number of Bloom filter bits per k-mer occurring at least twice in the reads (default is 14)
  -o, --output=STRING         Filename for output
  -c, --chunk-size=INT        Read chunksize to split betweeen threads (default is 10000)
  -v, --verbose               Print lots of messages during run
```

To obtain quickly the arguments *-n* and *-N* from a FASTQ file (if you have no idea), the tool KmerStream can be used (https://github.com/pmelsted/KmerStream). KmerStream output 3 numbers *F0*, *f1* and *F1*. You can then configure BFGraph with *F0* for *-n* and *F0-f1* for *-N*.
Arguments *-b* and *-B* basically control the Bloom filter false positive rates. A larger number means less false positives in the Bloom filters but more memory used. We advise to not modify these parameters unless you know what you are doing with Bloom filters.

### Building

```
./BFGraph contigs
```

```
Create a compacted de Bruijn graph from filtered FASTA/FASTQ files and save it to a GFA file

Usage: BFGraph contigs [arguments] ... FASTQ_files

Required arguments:
  -t, --threads=INT           Number of threads (default is 1)
  -k, --kmer-size=INT         Size of k-mers, same as for filtering reads (default is 31)
  -g, --min-size=INT          Size of minimizers, same as for filtering reads (default is 23)
  -f, --filtered=STRING       Filtered reads file
  -o, --output=STRING         Prefix for output files
  -c, --chunk-size=INT        Read chunksize to split betweeen threads (default is 10000)

Optional arguments:
  -n, --clip-tips             Clip tips shorter than k k-mers in length
  -d, --rm-isolated           Delete isolated contigs shorter than k k-mers in length
  -v, --verbose               Print lots of messages during run
```

## Notes

* BFGraph has been developed on x86-64 GNU/Linux and is compatible with MacOS. Porting it to a non-unix platform should be relatively easy, but we have not done it so far.

* If you run into bugs/problems or have questions/suggestions, please feel free to file an issue on this GitHub repository.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The kseq functions for reading fast(a|q)(.gz) files are copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

*   This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
