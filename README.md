# Bifrost

Highly Parallel and Memory Efficient Colored and Compacted de Bruijn Graph Construction

This repository contains the source code for a new parallel and memory efficient algorithm enabling the direct construction of the compacted de Bruijn graph without producing the intermediate uncompacted de Bruijn graph. Despite making extensive use of a probabilistic data structure (the Bloom filter), our algorithm guarantees that the produced compacted de Bruijn graph is deterministic. Furthermore, the algorithm features de Bruijn graph simplification steps used by assemblers such as tip clipping and isolated unitig removal. In addition, as disk-based software performance is significantly affected by the discrepancy of speed among disk storage technologies, our method uses only main memory storage.

## Dependencies

In order to compile and use Bifrost, you need a machine running a 64 bits POSIX-compliant UNIX or MacOS operating system. Bifrost successfully compiles and runs on Ubuntu 17.04 and MacOS.

In order to compile Bifrost, you will need:

- Cmake (https://cmake.org/)
- Jemalloc (http://www.canonware.com/jemalloc)
- Roaring (https://github.com/RoaringBitmap/CRoaring)
- zlib (https://zlib.net/)

All can be downloaded and installed by following the instructions on their respective websites. It is however most likely that at least few of them are available via a package manager for your operating system.

If you operating system is Ubuntu/Debian, you can install Cmake-Jemalloc-Zlib as follows:
```
sudo apt-get install cmake libjemalloc1 libjemalloc-dev zlib1g
```

If you operating system is MacOS, Cmake-Jemalloc-Zlib can be easily downloaded and installed via Homebrew:
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
cd <bifrost_directory>
mkdir build
cd build
cmake ..
make
make install
```

`make install` might requires `sudo` (`sudo make install`) to proceed.

By default, the installation creates a binary (*Bifrost*), a dynamic library (*libbifrost.so* for Unix or *libbifrost.dylib* for MacOS) and a static library (*libbifrost.a*).

The default maximum *k*-mer size supported is 31. To work with larger *k*, you must replace *MAX_KMER_SIZE* in *CMakeLists.txt* with a larger (appropriate) number, such as:
```
set( MAX_KMER_SIZE "64")
```

In this case, the maximum *k* allowed is 63. Keep in mind that increasing MAX_KMER_SIZE increases Bifrost memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

## Binary usage:

Type the following command to run Bifrost:
```
Bifrost
```

It should display the command line interface:
```
Bifrost x.y

Highly Parallel and Memory Efficient Colored and Compacted de Bruijn Graph Construction

Usage: Bifrost [Parameters] -o <output_prefix> -f <file_1> ...

Mandatory parameters with required argument:

  -f, --input-files        Input sequence files (FASTA or FASTQ, possibly gziped) and/or graph files (GFA)
  -o, --output-file        Prefix for output file (GFA output by default)

Optional parameters with required argument:

  -t, --threads            Number of threads (default is 1)
  -k, --kmer-length        Length of k-mers (default is 31)
  -g, --min-length         Length of minimizers (default is 23)
  -n, --num-kmers          Estimated number of different k-mers in input files (default: KmerStream estimation)
  -N, --num-kmers2         Estimated number of different k-mers occurring twice or more in the input files (default: KmerStream estimation)
  -b, --bloom-bits         Number of Bloom filter bits per k-mer occurring at least once in the input files (default is 14)
  -B, --bloom-bits2        Number of Bloom filter bits per k-mer occurring at least twice in the input files (default is 14)
  -l, --load-mbbf          Filename for input Blocked Bloom Filter, skips filtering step (default is no input)
  -w, --write-mbbf         Filename for output Blocked Bloom Filter (default is no output)
  -s, --chunk-size         Read chunksize to split between threads (default is 10000)

Optional parameters with no argument:

  -p, --produce-colors     Produce a colored and compacted de Bruijn graph
  -r, --reference          Reference mode, no filtering
  -i, --clip-tips          Clip tips shorter than k k-mers in length
  -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length
  -m, --keep-mercy         Keep low coverage k-mers connecting tips
  -a, --fasta              Output file is in FASTA format (only sequences) instead of GFA
  -v, --verbose            Print information messages during construction
```

Bifrost works in two steps:

1. reads are filtered to remove errors
2. the compacted de Bruijn graph is built from the filtered reads

If you want to input assembled genomes, use the `-r` parameter and no filtering will be applied, all k-mers of the files will be used to build the graph.

### Without colors

By default, Bifrost produces a compacted de Bruijn graph without colors.

### With colors (pre-alpha)

Colors are used to annotate k-mers with the set of genomes/samples in which they occur. Producing a colored and compacted de Bruijn graph using Bifrost is a two steps process:

1. **Build the compacted de Bruijn graph for each color (sample(s) or genome(s)) and output the unitigs of each graph to a FASTA file (`-a` parameter)**

   Example: Three input files *A.fastq*, *B.fastq* and *C.fastq*: A and B are reads to group in one color, C is an assembled genome to associate with another color.
   ```
   Bifrost -k 31 -t 4 -i -d -a -o AB_cdBG -f A.fastq B.fastq
   Bifrost -k 31 -t 4 -r -a -o C_cdBG -f C.fastq 
   ```
   In this example, each compacted de Bruijn graph is built using 31-mers (`-k 31`) and 4 threads (`-t 4`). For the read files A and B (`-f A.fastq B.fastq`), graph simplification steps are performed after construction (`-i -d`) and the graph is output to the FASTA file *AB_cdBG.fasta* (`-a -o AB_cdBG`). For the assembled genome C, the graph is built in reference mode (`-r`) and output to the FASTA file *C_cdBG.fasta* (`-a -o C_cdBG`). 

2. **Build the colored and compacted de Bruijn graph (`-p` parameter) using the previously produced FASTA files in input**

   Example: Two input FASTA files *AB_cdBG.fasta* and *C_cdBG.fasta*. Each file contains the unitigs of a compacted de Bruijn graph and is going to be represented by a color.
   ```
   Bifrost -k 31 -t 4 -p -o ABC_ccdBG -f AB_cdBG.fasta C_cdBG.fasta
   ```
   In this example, the colored (`-p`) and compacted de Bruijn graph is built using 31-mers (`-k 31`) and 4 threads (`-t 4`) from the files *AB_cdBG.fasta* and *C_cdBG.fasta* (`-f AB_cdBG.fasta C_cdBG.fasta`). The graph will be output to a GFA file *ABC_ccdBG.gfa* and colors will be output to file *ABC_ccdBG.bfg_colors* (`-o ABC_ccdBG`).

<img src="pipeline_colored_cdbg.png" alt="pipeline_colors" style="width: 200px;"/>

02-02-2018: More color options coming soon

## API

(Work in progress)

### Documentation

Documentation for the Bifrost library is available in the */doc/doxygen* folder (HTML version, open *index.html*).

The following command regenerates the documentation:
```
cd <bifrost_directory>
doxygen Doxyfile
```

The documentation contains a description of all the functions and structures of the library.

TODO: Code snippets

### Usage

Once Bifrost is installed on your operating system, just use
```
#include <bifrost/CompactedDBG.hpp>
```
in your C++ code. Then, use the following flags for linking:
```
-lbifrost -pthread
```

`-lbifrost` refers to the Bifrost dynamic library which is multi-threaded and hence, requires a threading library, usually the POSIX Threads library `-pthread`. It is possible while compiling your program with the Bifrost library that your compiler complains about missing dependencies. In that case, use the following flags for linking:

```
-lbifrost -ljemalloc -lroaring -pthread -lz
```

You can also link to the Bifrost static library (*libbifrost.a*) for better performance:
```
<path_to_lib_folder>/libbifrost.a -ljemalloc -lroaring -pthread -lz
```

### With colors (pre-alpha)

```
#include <bifrost/ColoredCDBG.hpp>
```

02-02-2018: No documentation available yet. As a pre-alpha version, the API might be changing in the future and it might still contain some bugs.

## FAQ

**What are the accepted input file formats?**

FASTA , FASTQ and GFA. Input FASTA and FASTQ files can be compressed with gzip (extension .gz). If you input a GFA file, you probably want to run Bifrost in reference mode (`-r` parameter, build the graph from all k-mers of the sequences).

**Can I mix different file formats in input?**

Yes, as long as they are FASTA, FASTQ and GFA.

**If I input a GFA file, does it need to contain already a compacted de Bruijn graph?**

No, it can contain any type of sequence graph (like an uncompacted de Bruijn graph). Bifrost will extract the sequences from the file and build the compacted de Bruijn graph out of them.

**Can I build a compacted de Bruijn graph from assembled genomes and reads?**

Yes. First, run Bifrost with your assembled genomes only (`-r` parameter) and output the unitigs to a FASTA file (`-a` parameter). Then, run Bifrost a second time with the previously produced FASTA file and your read files.

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

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
