# Bifrost

### Highly parallel construction and indexing of colored and compacted de Bruijn graphs

* **Build**, **index** and **color** the compacted de Bruijn graph
* **No need to build the uncompacted** de Bruijn graph
* **Reads** or **assembled genomes** as input
* Output **graph in GFA** (can be visualized with [Bandage](https://github.com/rrwick/Bandage))
* **Graph cleaning**: short tip clipping, etc.
* **No disk** usage (adapted for cluster architectures)
* **Multi-threaded** and **SIMD** optimized
* **No parameters to estimate** with other tools
* **C++ API** available:
    * Associate **your data with vertices**
    * **Add** or **remove** (sub-)sequences / *k*-mers / colors
    * **Find unitigs** containing **queried k-mers**

## Requirements

* 64 bits POSIX-compliant UNIX or MacOS operating system
* C++11 capable compiler:
    * [GCC](https://gcc.gnu.org/) 5.0 or later
    * [Clang](http://clang.llvm.org/) 3.5 or later
* [Cmake](https://cmake.org/) 2.8.12 or later
* [Zlib](https://zlib.net/)

GCC/Clang, Cmake and Zlib are probably already installed on your computer (those are installed by default on most operating systems) but they can be downloaded and installed by following the instructions on their respective websites. It is most likely that all are available via a package manager for your operating system: 

* Ubuntu/Debian:
```
sudo apt-get install build-essential cmake zlib1g
```
* MacOS (with [Homebrew](https://brew.sh/)):
```
brew install --with-toolchain llvm
brew install cmake zlib
```

## Compilation and Installation

```
cd <bifrost_directory>
mkdir build
cd build
cmake ..
make
make install
```

`make install` might requires `sudo` (`sudo make install`) to proceed. See [Troubleshooting](#troubleshooting) if you have any problem during the installation.

By default, the installation creates:
* a binary (*Bifrost*)
* a dynamic library (*libbifrost.so* for Unix or *libbifrost.dylib* for MacOS)
* a static library (*libbifrost.a*).

The default maximum *k*-mer size supported is 31. To work with larger *k*, you must replace *MAX_KMER_SIZE* in *CMakeLists.txt* with a larger (appropriate) number, such as:
```
set( MAX_KMER_SIZE "64")
```

In this case, the maximum *k* allowed is 63. Keep in mind that increasing *MAX_KMER_SIZE* increases Bifrost memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

## Binary usage:

```
Bifrost
```

displays the command line interface:
```
Bifrost x.y

Highly parallel construction and indexing of colored and compacted de Bruijn graphs

Usage: Bifrost [COMMAND] [GENERAL_PARAMETERS] [COMMAND_PARAMETERS]

[COMMAND]:

   build                   Build a compacted de Bruijn graph, with or without colors
   update                  Update a compacted (possible colored) de Bruijn graph with new sequences

[GENERAL_PARAMETERS]:

   > Mandatory with required argument:

   -s, --input-seq-files    Input sequence files (FASTA/FASTQ possibly gzipped)
                            Input sequence files can be provided as a list in a TXT file (one file per line)
                            K-mers with exactly 1 occurrence in the input sequence files will be discarded
   -r, --input-ref-files    Input reference files (FASTA/FASTQ possibly gzipped and GFA)
                            Input reference files can be provided as a list in a TXT file (one file per line)
                            All k-mers of the input reference files are used
   -o, --output-file        Prefix for output file(s)

   > Optional with required argument:

   -t, --threads            Number of threads (default is 1)

   > Optional with no argument:

   -i, --clip-tips          Clip tips shorter than k k-mers in length
   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length
   -v, --verbose            Print information messages during execution

[COMMAND_PARAMETERS]: build

   > Optional with required argument:

   -k, --kmer-length        Length of k-mers (default is 31)
   -m, --min-length         Length of minimizers (default is 23)
   -b, --bloom-bits         Number of Bloom filter bits per k-mer with 1+ occurrences in the input files (default is 14)
   -B, --bloom-bits2        Number of Bloom filter bits per k-mer with 2+ occurrences in the input files (default is 14)
   -l, --load-mbbf          Input Blocked Bloom Filter file, skips filtering step (default is no input)
   -w, --write-mbbf         Output Blocked Bloom Filter file (default is no output)
   -u, --chunk-size         Read chunk size per thread (default is 64)

   > Optional with no argument:

   -c, --colors             Color the compacted de Bruijn graph (default is no coloring)
   -y, --keep-mercy         Keep low coverage k-mers connecting tips
   -a, --fasta              Output file is in FASTA format (only sequences) instead of GFA

[COMMAND_PARAMETERS]: update

   > Mandatory with required argument:

   -g, --input-graph-file   Input graph file to update (GFA format)

   > Optional with required argument:

   -f, --input-color-file   Input color file associated with the input graph file to update
   -k, --kmer-length        Length of k-mers (default is read from input graph file if built with Bifrost or 31)
   -m, --min-length         Length of minimizers (default is read from input graph file if built with Bifrost or 23)
```

### Examples

1. **Build a compacted de Bruijn graph from read files and clean the graph**
   ```
   Bifrost build -t 4 -k 31 -i -d -o AB_graph -s A.fastq B.fastq
   ```
   The compacted de Bruijn graph is built (`build`) with 4 threads (`-t 4`) from the 31-mers (`-k 31`) of files *A.fastq* and *B.fastq* (`-s A.fastq B.fastq`). By using parameter `-s`, files *A.fastq* and *B.fastq* are filtered: 31-mers occurring exactly once in *A* and *B* are discarded from the construction. Graph simplification steps are performed after building (`-i -d`) and the graph is written to file *AB_graph.gfa* (`-o AB_graph`).

2. **Build a compacted de Bruijn graph from a reference genome file**
   ```
   Bifrost build -t 4 -k 31 -o C_graph -r C.fasta
   ```
   The compacted de Bruijn graph is built (`build`) with 4 threads (`-t 4`) from the 31-mers (`-k 31`) of file *C.fasta* (`-r C.fasta`). By using parameter `-r`, file *C.fasta* is NOT filtered: all 31-mers occurring in *C* are used during the construction. The graph is written to file *C_graph.gfa* (`-o C_graph`).

3. **Building a compacted and colored de Bruijn graph from read files and reference genome files, clean the graph**
   ```
   Bifrost build -t 4 -k 31 -c -i -d -o ABC -s A.fastq B.fastq -r C.fasta
   ```
   Combining the two previous examples, the compacted de Bruijn graph is built (`build`) with 4 threads (`-t 4`) from the 31-mers (`-k 31`) of files *A.fastq*, *B.fastq* (`-s A.fastq B.fastq`) and file *C.fasta* (`-r C.fasta`). Graph simplification steps are performed after building (`-i -d`). The graph is colored (`-c`), meaning that each k-mer of the graph unitigs keeps track of whether it occurs in *A*, *B* or *C*. The graph is written to file *ABC.gfa* and the colors are written to file *ABC.bfg_colors* (`-o ABC`).

4. **Updating a compacted de Bruijn graph with a reference genome file**
   ```
   Bifrost update -t 4 -o CD_graph -r D.fasta -g C_graph.gfa
   ```
   The compacted de Bruijn graph *C* (`-g C_graph.gfa`) is updated (`update`) with 4 threads (`-t 4`) from the *k*-mers of file *D.fasta* (`-r D.fasta`). By using parameter `-r`, file *D.fasta* is NOT filtered: all *k*-mers occurring in *D* are used during the merging. The graph is written to file *CD_graph.gfa* (`-o CD_graph`).

5. **Updating a compacted and colored de Bruijn graph with read files and clean the graph**
   ```
   Bifrost update -t 4 -i -d -o ABCEF -s E.fastq F.fastq -g ABC.gfa -f ABC.bfg_colors
   ```
   The compacted and colored de Bruijn graph *ABC* (`-g ABC.gfa -f ABC.bfg_colors`) is updated (`update`) with 4 threads (`-t 4`) from the *k*-mers of files *E.fastq* and *F.fastq* (`-s E.fastq F.fastq`). Graph simplification steps are performed after merging (`-i -d`). The graph is written to file *ABCEF.gfa* and the colors are written to file *ABCEF.bfg_colors* (`-o ABCEF`).

## API

(Work in progress)

### Documentation

Documentation for the Bifrost library is available in the */doc/doxygen* folder (HTML version, open *html/index.html*).

The following command regenerates the documentation:
```
cd <bifrost_directory>
doxygen Doxyfile
```

The documentation contains a description of all the functions and structures of the library.

### Usage

Once Bifrost is installed on your operating system, just use
```
#include <bifrost/CompactedDBG.hpp>
```
in your C++ code. Then, use the following flags for compiling:
```
-O3 -std=c++11 -march=native
```

and the following flags for linking:
```
-lbifrost -pthread -lz
```

You can also link to the Bifrost static library (*libbifrost.a*) for better performance:
```
<path_to_lib_folder>/libbifrost.a -pthread -lz
```

### With colors (beta)

```
#include <bifrost/ColoredCDBG.hpp>
```

### Changelog

* **08-07-2018**
	* Add de Bruijn graphs merging functions (`CompactedDBG::merge()` and `ColoredCDBG::merge()`) and addition assignment operators (`CompactedDBG::operator+=()` and `ColoredCDBG::operator+=()`, same as `merge()` but uses only one thread).
	* Add de Bruijn graphs comparison functions `CompactedDBG::operator==()`, `CompactedDBG::operator!=()`, `ColoredCDBG::operator==()` and `ColoredCDBG::operator!=()`.
	* Delete `CompactedDBG::empty()` and `ColoredCDBG::empty()` to be consistent with STD containers (those functions were emptying the graph of its content while `empty()` of STD containers returns whether the container is empty). Now, to empty the graph, use `CompactedDBG::clear()` and `ColoredCDBG::clear()`.
    * Major changes in the abstract class `CDBG_Data_t` and `CCDBG_Data_t`:
    	* All the functions are now **non**-static.
    	* Function `join()` is renamed `concat()` and works a bit differently (have a look at the doc). Quickly, `join()` was concatenating two unitigs A and B such that the result was A=AB and B was deleted from the graph. Now, `concat()` deletes A and B from the graph and adds a new unitig C=AB.
    	* Function `sub()` is renamed `extract()`.
    	* Add the functions `merge()` and `clear()` which **must** be overloaded too in the derived class of `CDBG_Data_t` and `CCDBG_Data_t`.
    * Bugfix with graph reading from disk and concatenating UnitigColors.

## FAQ

**What are the accepted input file formats?**

FASTA, FASTQ and GFA. Input FASTA and FASTQ files can be compressed with gzip (extension .gz). If you input a GFA file for the construction, you probably want to use the `-r` parameter for that file.

**Can I mix different file formats in input?**

Yes, as long as they are FASTA, FASTQ and GFA.

**Can I provide in input a file which is a list of files?**

Yes, a text file (which is a list of files) can be given in input. It must contain one input filename per line with no empty lines.

**If I input a GFA file, does it need to contain an already compacted de Bruijn graph?**

No, it can contain any type of sequence graph (like an uncompacted de Bruijn graph or a sequence graph).

**Can I build a compacted (colored) de Bruijn graph from assembled genomes and reads?**

Yes. Input your assembled genomes with parameter `-r` and your reads with parameter `-s`.

**In which order are inserted the colors?**

A color corresponds to an input file the graph was built/updated from. The order in which the colors are inserted is the same as the order of the files given by parameter `-r` and parameter `-s`. However, in case both parameters `-r` and `-s` are used, no assumption can be made on whether the files given by parameter `-s` will be inserted before or after the ones given by parameter `-r`.

## Troubleshooting

The following might happen when environment variables are not set correctly on your system:

* compilation (`make`) fails because some header files (*.h*) are not found

Assuming the header files (*.h*) are located at the path */usr/local/include/*, the following command set the environment variables *C_INCLUDE_PATH* and *CPLUS_INCLUDE_PATH* correctly for the time of the session:
```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include/
```

* executing the binary of Bifrost fails because *libbifrost.so* or *libbifrost.a* is not found

Assuming that *libbifrost*.(*so*|*dylib*|*a*) is located at the path */usr/local/lib/*, the following command set the environment variables *LD_LIBRARY_PATH*, *LIBRARY_PATH* and *PATH* correctly for the time of the session:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/lib/
export PATH=$PATH:/usr/local/lib/
```

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.

## License

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

* The kseq library is copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)

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
