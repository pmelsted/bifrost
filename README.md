# Bifrost

### Parallel construction, indexing and querying of colored and compacted de Bruijn graphs

* **Build**, **index**, **color** and **query** the compacted de Bruijn graph
* **No need to build the uncompacted** de Bruijn graph
* **Reads** or **assembled genomes** as input
* Output **graph in GFA** (can be visualized with [Bandage](https://github.com/rrwick/Bandage)), **FASTA** or **binary**
* **Graph cleaning**: short tip clipping, etc.
* **Multi-threaded**
* **No parameters to estimate** with other tools
* **Exact** or **approximate** *k*-mer search of queries
* **C++ API** available:
    * Associate **your data with vertices**
    * **Add** or **remove** (sub-)sequences / *k*-mers / colors
    * **Find unitigs** containing **queried k-mers**

## Table of Contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Binary usage](#binary-usage)
* [API](#api)
* [FAQ](#faq)
* [Troubleshooting](#troubleshooting)
* [Citation](#citation)
* [Contact](#contact)
* [License](#license)

## Requirements

To install Bifrost using Bioconda, go directly to Section [Installation](#installation). To install from source, you will need:

* C++11 compiler:
    * [GCC](https://gcc.gnu.org/) >= 5.1.0
    * [Clang](http://clang.llvm.org/) >= 3.5
* [Cmake](https://cmake.org/) >= 2.8.12
* [Zlib](https://zlib.net/)

All are probably already installed on your computer as those are installed by default on most operating systems. They can be downloaded and installed by following the instructions on their respective websites. However, it is most likely they are all available via a package manager for your operating system: 

* **Ubuntu/Debian**:
```
sudo apt-get install build-essential cmake zlib1g-dev
```
* **MacOS** (with [Homebrew](https://brew.sh/)):
```
brew install --with-toolchain llvm
brew install cmake zlib
```
* **Windows 10**:

1. Open the Windows Store
2. Search and install the `Ubuntu` app (from `Canonical Group Limited`)
3. Open the Windows main menu and open the `Ubuntu` app (it should open an Ubuntu terminal)
4. Use the following command in the Ubuntu terminal:
```
sudo apt-get install build-essential cmake zlib1g-dev
```
5. Use the opened Ubuntu terminal for compiling, installing and running Bifrost (see next section). See [Troubleshooting](#troubleshooting) if you have any problem during the installation.

## Installation

Compared to the source install, the Conda package might not have the latest Bifrost version, does not support *k>31* nor native compilation. Use the source installation for benchmarking.

* From [Bioconda](https://bioconda.github.io):

  ```
  conda install -c bioconda bifrost
  ```

* From source

  ```
  git clone https://github.com/pmelsted/bifrost.git
  cd bifrost && mkdir build && cd build
  cmake ..
  make
  make install
  ```

  `make install` might require `sudo` (`sudo make install`) to proceed. To install Bifrost in the non-default path `/some/path/`, add the option `-DCMAKE_INSTALL_PREFIX=/some/path/` to the `cmake` command.

  By default, the installation creates:
  * a binary (*Bifrost*)
  * a dynamic library (*libbifrost.so* for Unix or *libbifrost.dylib* for MacOS)
  * a static library (*libbifrost.a*)

  **Advanced options**
  * Bifrost compiles by default with `-march=native`: the compiler targets architecture instructions specific to the machine Bifrost is compiled on. Hence, the binary and library produced might not work on a different machine. Native compilation can be disabled by adding the option `-DCOMPILATION_ARCH=OFF` to the `cmake` command (disables all AVX2 optimizations too). Alternatively, you can use this option to specify the architecture you want to target: `x86-64`, `knl`, etc. Default is `-DCOMPILATION_ARCH=native`.
  * Bifrost uses AVX2 instructions during graph construction which can be disabled by adding the option `-DENABLE_AVX2=OFF` to the `cmake` command.

  If you encounter any problem during the installation, see the [Troubleshooting](#troubleshooting) section.

### Large *k*-mers

The default maximum *k*-mer size supported is 31. To work with larger *k* in the binary, you must install Bifrost from source and replace *MAX_KMER_SIZE* with a larger multiple of 32. This can be done in two ways:

* By adding the following option to the `cmake` command:
```
-DMAX_KMER_SIZE=64
```

* By replacing *MAX_KMER_SIZE* in *CMakeLists.txt*:
```
SET(MAX_KMER_SIZE "64" CACHE STRING "MAX_KMER_SIZE")
```

Actual maximum k-mer size is *MAX_KMER_SIZE-1*, e.g maximum *k* is 63 for *MAX_KMER_SIZE=64*. Increasing *MAX_KMER_SIZE* increases Bifrost memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

The maximum size of minimizers (*g*-mers) *MAX_GMER_SIZE* can be adjusted the same way as *MAX_KMER_SIZE*. This is especially useful if you want to use a large *k*-mer size but a small *g*-mer size. By default, *MAX_GMER_SIZE* is equal to *MAX_KMER_SIZE*.

To work with larger *k* when using the Bifrost API, the new value *MAX_KMER_SIZE* must be given to the compiler and linker as explained in Section [API](#api)

## Binary usage:

```
Bifrost
```

displays the command line interface:
```
Bifrost x.y.z

Highly parallel construction, indexing and querying of colored and compacted de Bruijn graphs

Usage: Bifrost [COMMAND] [PARAMETERS]

[COMMAND]:

   build                   Build a compacted de Bruijn graph, with or without colors
   update                  Update a compacted (colored) de Bruijn graph with new sequences
   query                   Query a compacted (colored) de Bruijn graph

[PARAMETERS]: build

   > Mandatory with required argument:

   -s, --input-seq-file     Input sequence file in fasta/fastq(.gz) format
                            Multiple files can be provided as a list in a text file (one file per line)
                            K-mers with exactly 1 occurrence in the input sequence files will be discarded
   -r, --input-ref-file     Input reference file in fasta/fastq(.gz) or gfa(.gz) format
                            Multiple files can be provided as a list in a text file (one file per line)
                            All k-mers of the input reference files are used
   -o, --output-file        Prefix for output file(s)

   > Optional with required argument:

   -t, --threads            Number of threads (default: 1)
   -k, --kmer-length        Length of k-mers (default: 31)
   -m, --min-length         Length of minimizers (default: auto)
   -B, --bloom-bits         Number of Bloom filter bits per k-mer (default: 24)
   -T, --tmp-dir            Path for tmp directory (default: creates tmp directory in output directory)
   -l, --load-mbbf          Input Blocked Bloom Filter file, skips filtering step (default: no input)
   -w, --write-mbbf         Output Blocked Bloom Filter file (default: no output)

   > Optional with no argument:

   -c, --colors             Color the compacted de Bruijn graph
   -i, --clip-tips          Clip tips shorter than k k-mers in length
   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length
   -f, --fasta-out          Output file in fasta format (only sequences) instead of gfa (unless graph is colored)
   -b, --bfg-out            Output file in bfg/bfi format (Bifrost graph/index) instead of gfa (unless graph is colored)
   -n, --no-compress-out    Output files must be uncompressed
   -N, --no-index-out       Do not make index file
   -v, --verbose            Print information messages during execution

[PARAMETERS]: update

  > Mandatory with required argument:

   -g, --input-graph-file   Input graph file to update in gfa(.gz) or bfg format
   -s, --input-seq-file     Input sequence file in fasta/fastq(.gz) format
                            Multiple files can be provided as a list in a text file (one file per line)
                            K-mers with exactly 1 occurrence in the input sequence files will be discarded
   -r, --input-ref-file     Input reference file in fasta/fastq(.gz) or gfa(.gz) format
                            Multiple files can be provided as a list in a text file (one file per line)
                            All k-mers of the input reference files are used
   -o, --output-file        Prefix for output file(s)

   > Optional with required argument:

   -I, --input-index-file   Input index file associated with graph to update in bfi format
   -C, --input-color-file   Input color file associated with graph to update in color.bfg format
   -t, --threads            Number of threads (default: 1)
   -k, --kmer-length        Length of k-mers (default: read from input graph file if built with Bifrost or 31)
   -m, --min-length         Length of minimizers (default: read from input graph if built with Bifrost, auto otherwise)
   -T, --tmp-dir            Path for tmp directory (default: creates tmp directory in output directory)

   > Optional with no argument:

   -i, --clip-tips          Clip tips shorter than k k-mers in length
   -d, --del-isolated       Delete isolated contigs shorter than k k-mers in length
   -f, --fasta-out          Output file in fasta format (only sequences) instead of gfa (unless colors are output)
   -b, --bfg-out            Output file in bfg/bfi format (Bifrost graph/index) instead of gfa (unless graph is colored)
   -n, --no-compress-out    Output files must be uncompressed
   -N, --no-index-out       Do not make index file
   -v, --verbose            Print information messages during execution

[PARAMETERS]: query

  > Mandatory with required argument:

   -g, --input-graph-file   Input graph file to query in gfa(.gz) or bfg
   -q, --input-query-file   Input query file in fasta/fastq(.gz). Each record is a query.
                            Multiple files can be provided as a list in a text file (one file per line)
   -o, --output-file        Prefix for output file

   > Optional with required argument:

   -e, --min_ratio-kmers    Minimum ratio of k-mers from each query that must occur in the graph
   -E, --min-nb-colors      Minimum number of colors from each query that must occur in the graph
   -I, --input-index-file   Input index file associated with graph to query in bfi format
   -C, --input-color-file   Input color file associated with the graph to query in color.bfg format
   -t, --threads            Number of threads (default: 1)
   -k, --kmer-length        Length of k-mers (default: read from input graph if built with Bifrost or 31)
   -m, --min-length         Length of minimizers (default: read from input graph if built with Bifrost, auto otherwise)
   -T, --tmp-dir            Path for tmp directory (default: creates tmp directory in output directory)

   > Optional with no argument:

   -Q, --files-as-queries   All fastq/fastq records in each input query file constitute a single query.
   -p, --ratio-found-km     Output the ratio of found k-mers for each query (disable parameters -e and -E)
   -a, --approximate        Graph is searched using exact and inexact k-mers (1 substitution or indel allowed per k-mer)
   -v, --verbose            Print information messages during execution
```

### Use cases

The following use cases describe some simple and common usage of the Bifrost CLI. However, many options are provided by the CLI to perform more specific actions (graph cleaning, approximate querying, etc.).

- **Build**

  1. **Build a compacted de Bruijn graph from read files**
     ```
     Bifrost build -t 4 -k 31 -s A.fastq -s B.fastq -o AB
     ```
     The compacted de Bruijn graph is built (`build`) using 4 threads (`-t 4`) from the 31-mers (`-k 31`) of files *A.fastq* and *B.fastq* (`-s A.fastq -s B.fastq`). By using parameter `-s`, files *A.fastq* and *B.fastq* are filtered: 31-mers occurring exactly once in *A* and *B* are discarded from the construction. The graph is written to file *AB.gfa.gz* and a Bifrost index is written to file *AB.bfi* (`-o AB`).

  2. **Build a compacted de Bruijn graph from a reference genome file**
     ```
     Bifrost build -t 4 -k 31 -r C.fasta -o C
     ```
     Same as previous use case but by using parameter `-r`, file *C.fasta* is NOT filtered: all 31-mers occurring in *C* are used during the construction. The graph is written to file *C.gfa.gz* and a Bifrost index is written to file *C.bfi* (`-o C`).

  3. **Build a compacted and colored de Bruijn graph from read files and reference genome files**
     ```
     Bifrost build -t 4 -k 31 -c -s A.fastq -s B.fastq -r C.fasta -o ABC 
     ```
     Combining the two previous use cases, the compacted de Bruijn graph is built (`build`) using 4 threads (`-t 4`) from the filtered 31-mers (`-k 31`) of files *A.fastq* and *B.fastq* (`-s A.fastq -s B.fastq`) and the unfiltered 31-mers of file *C.fasta* (`-r C.fasta`). The graph is colored (`-c`) such that for each k-mer in the unitigs of the graph is recorded whether it occurs in *A*, *B* or *C*. The graph is written to file *ABC.gfa.gz*, its colors are written to file *ABC.color.bfg* and a Bifrost index is written to file *ABC.bfi* (`-o ABC`).

     Additional options of interest for building are:
     - `-i`: Delete all tips composed of unitigs shorter than *k* *k*-mers
     - `-d`: Delete all connected components composed of one unitig shorter than *k* *k*-mers

- **Update**

  1. **Update a compacted de Bruijn graph with a reference genome file**
     ```
     Bifrost update -t 4 -g A_graph.gfa.gz -r B.fasta -o AB
     ```
     The compacted de Bruijn graph *A* (`-g A_graph.gfa.gz`) is updated (`update`) using 4 threads (`-t 4`) with the unfiltered *k*-mers of file *B.fasta* (`-r B.fasta`). The Bifrost index *A_graph.bfi* is automatically loaded if available in the same path as the graph but can also be loaded with `-I`. The graph is written to file *AB.gfa.gz* and a Bifrost index is written to file *AB.bfi* (`-o AB`).

  2. **Update a compacted and colored de Bruijn graph with read files**
     ```
     Bifrost update -t 4 -g A.gfa.gz -f A.color.bfg -s B.fastq -s C.fastq -o ABC
     ```
     The compacted and colored de Bruijn graph *A* (`-g A.gfa.gz -f A.color.bfg`) is updated (`update`) using 4 threads (`-t 4`) with the filtered *k*-mers of files *B.fastq* and *C.fastq* (`-s B.fastq -s C.fastq`). The Bifrost index *A.bfi* is automatically loaded if available in the same path as the graph but can also be loaded with `-I`. The merged graph is written to file *ABC.gfa.gz*, its colors are written to file *ABC.color.bfg* and a Bifrost index is written to file *ABC.bfi* (`-o ABC`).

        Additional options of interest for merging are:
     - `-i`: Delete all tips composed of unitigs shorter than *k* *k*-mers
     - `-d`: Delete all connected components composed of one unitig shorter than *k* *k*-mers

- **Query**

  The default querying behavior is to report the number of *k*-mers shared between the queries and:
  - the graph if input graph is **not** colored
  - each color if input graph is colored
 
  1. **Query a compacted de Bruijn graph for the number of k-mers shared between the queries and the graph**
     ```
     Bifrost query -t 4 -g A.gfa.gz -q in_queries.fasta -o out_queries_result 
     ```
     The compacted de Bruijn graph *A* (`-g A.gfa.gz`) is queried (`query`) using 4 threads (`-t 4`) for the number of *k*-mers shared between the sequences of file *in_queries.fasta* (`-q in_queries.fasta`) and the graph *A*. The Bifrost index *A.bfi* is automatically loaded if available in the same path as the graph but can also be loaded with `-I`. The results are stored in a matrix written to file *out_queries_result.tsv* (`-o out_queries_result`): rows are the queries, column is the graph, intersection of row/column is the number of shared *k*-mers between the query and the graph.

  2. **Query a compacted de Bruijn graph for presence/absence of queries in the graph**
     ```
     Bifrost query -t 4 -e 0.8 -g A.gfa.gz -q in_queries.fasta -o out_queries_result 
     ```
     Same as previous use case but instead of returning a number of shared *k*-mers per query, it returns a binary value indicating whether the query is present in the graph or not (1 if present, 0 if absent). At least 80% of the *k*-mers in each query must be found in the graph to report the query as present (`-e 0.8`).

  3. **Query a colored and compacted de Bruijn graph for the number of k-mers shared between the queries and the colors of the graph**
     ```
     Bifrost query -t 4 -g AB.gfa.gz -C AB.color.bfg -q in_queries.fasta -o out_queries_result 
     ```
     The compacted and colored de Bruijn graph *AB* (`-g AB.gfa.gz -C AB.color.bfg`) is queried (`query`) using 4 threads (`-t 4`) for the number of *k*-mers shared between the sequences of file (`-q in_queries.fasta`). The Bifrost index *AB.bfi* is automatically loaded if available in the same path as the graph but can also be loaded with `-I`. The results are stored in a matrix written to file *out_queries_result.tsv* (`-o out_queries_result`): rows are the queries, columns are the colors, intersection of row/column is an integer indicating the number of *k*-mers from the query occuring in the graph with the corresponding color.

  4. **Query a colored and compacted de Bruijn graph for presence/absence of queries in each color of the graph**
     ```
     Bifrost query -t 4 -e 0.8 -g AB.gfa.gz -C AB.color.bfg -q in_queries.fasta -o out_queries_result 
     ```
     Same as previous use case but instead of returning a number of shared *k*-mers per query and color, it returns a binary value indicating whether the query is present in the graph with the corresponding color or not (1 if present, 0 if absent). At least 80% of the *k*-mers in each query must be found in the graph with corresponding color to report the query as present for that color (`-e 0.8`).

  Additional options of interest for querying are:
  -  `-p`: Outputs a ratio of shared *k*-mers (w.r.t the number of *k*-mers in each query) instead of the number of *k*-mers
  -  `-Q`: Performs the querying per file (with multiple sequences) rather than per sequence
  -  `-a`: Enable approximate *k*-mer matches
  -  `-E`: A *k*-mer is only reported as present if it is colored by *x* many colors
      
## API

Changes in the API are reported in the [Changelog](https://github.com/pmelsted/bifrost/blob/master/Changelog.md).

### Tutorial

The [API tutorial](doc/tutorial/Intro.md) should help you get started with the C++ API.

### Documentation

Documentation for the Bifrost library is available in the */doc/doxygen* folder (HTML version, open *html/index.html*).

The following command regenerates the documentation:
```
cd <bifrost_directory>
doxygen Doxyfile
```

The documentation contains a description of all the functions and structures of the library.

### Usage

The Bifrost C++ API can be used by adding
```
#include <bifrost/CompactedDBG.hpp>
```
for uncolored compacted de Bruijn graphs and
```
#include <bifrost/ColoredCDBG.hpp>
```
for colored compacted de Bruijn graphs in your C++ headers.

To compile, we recommend using the following compile flags:
```
-O3 -std=c++11
```
Furthermore, Bifrost compiles by default with flag `-march=native` so unless native compilation was disabled when installing Bifrost, use flag `-march=native` too.

Finally, use the following flags for linking:
```
-lbifrost -pthread -lz
```

You can also link to the Bifrost static library (*libbifrost.a*) for better performance:
```
<path_to_lib_folder>/libbifrost.a -pthread -lz
```

The default maximum *k*-mer size supported is 31. To work with larger *k*, the code using the Bifrost C++ API must be compiled and linked with the flag `-DMAX_KMER_SIZE=x` for compiling and linking where `x` is a larger multiple of 32, such as:
```
-DMAX_KMER_SIZE=64
```
Actual maximum k-mer size is *MAX_KMER_SIZE-1*, e.g maximum *k* is 63 for *MAX_KMER_SIZE=64*. Increasing *MAX_KMER_SIZE* increases Bifrost memory usage (*k*=31 uses 8 bytes of memory per *k*-mer while *k*=63 uses 16 bytes of memory per *k*-mer).

## FAQ

**Can I provide in input multiple files?**

Yes, use parameter `-r` or `-s` for each file to input.

**Can I provide in input a file which is a list of files?**

Yes, a text file containing one input filename per line with no empty lines can be given in input.

**What are the accepted input file formats?**

FASTA, FASTQ and GFA. Input FASTA and FASTQ files can be compressed with gzip (extension .gz). If you input a GFA file for the construction, you probably want to use the `-r` parameter for that file.

**Can I mix different file formats in input?**

Yes, as long as they are FASTA, FASTQ and GFA.

**If I input a GFA file for building the de Bruijn graph, does it need to contain an already compacted de Bruijn graph?**

No, it can contain any type of sequence graph (like an uncompacted de Bruijn graph or a sequence graph).

**Can I build a compacted (colored) de Bruijn graph from assembled genomes and reads?**

Yes. Input your assembled genomes with parameter `-r` and your reads with parameter `-s`.

**Can I use the graph file without its color file ?**

Yes. Just do not input the color file and Bifrost will consider it is an **un**colored compacted de Bruijn graph.

**In which order are inserted the colors?**

A color corresponds to an input file the graph was built/updated from. The order in which the colors are inserted is the same as the order of the files given by parameter `-r` and parameter `-s`. However, in case both parameters `-r` and `-s` are used, no assumption can be made on whether the files given by parameter `-s` will be inserted before or after the ones given by parameter `-r`.

**Different runs of Bifrost on the same dataset with the same parameters produces graphs with different unitigs. Which graph is correct?**

All of them. The difference between the graphs resides in circular unitigs (unitigs connecting to themselves) which are their own connected components ("isolated"). These unitigs can have a different sequence from one run to another because the starting position will be different, yet they represent exactly the same sequence. As an example, circular unitig ATAT composed of 3-mers can also be represented with sequence TATA. The number of unitigs will remain the same from one graph to another.

**Is it possible to get the colors per *k*-mer in a parsable (non-binary) file format?**

Yes, please see [this solution](https://github.com/pmelsted/bifrost/issues/50)

## Troubleshooting

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

* Bifrost crashes right at the beginning with error `Illegal instruction`

  You are most likely running Bifrost on a different machine than the one used to compile it. By default, Bifrost is compiled in native mode such to target architecture instructions specific to the machine it is compiled on. Using Bifrost on a different machine with a different architecture might result in this error. To solve this issue, Bifrost must be recompiled with native architecture compilation disabled, as explained in the Advanced options of Section [Installation](#installation).

## Citation

```
@article {holley2019bifrost,
  author = {Holley, Guillaume and Melsted, P{\'a}ll},
  title = "{Bifrost - Highly parallel construction and indexing of colored and compacted de Bruijn graphs}",
  elocation-id = {695338},
  doi = {10.1101/695338},
  journal = {bioRxiv},
  year = {2019}
}
```

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.

## License

* [Bifrost](https://github.com/pmelsted/bifrost/blob/master/LICENSE) is BSD2 licensed
* The [wyhash](https://github.com/wangyi-fudan/wyhash) library is Unlicense licensed
* The [popcount](https://github.com/kimwalisch/libpopcnt) library is BSD licensed
* The [fastmod](https://github.com/lemire/fastmod) library is Apache 2.0 licensed
* The [kseq](http://lh3lh3.users.sourceforge.net/kseq.shtml) library is copyrighted by Heng Li and released under the MIT license
* The [CRoaring](https://github.com/RoaringBitmap/CRoaring) library is Apache 2.0 licensed
* The [zstr](https://github.com/mateidavid/zstr) library is MIT licensed
* The GetRSS library is Creative Commons Attribution 3.0 licensed
