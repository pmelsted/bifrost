BFGraph is a memory efficient program for counting k-mers from sequencing files.

Installation
============
Just run make to compile.

% make

The default maximum k-mer size supported is 31 (8 bytes of memory per k-mer), to modify this either
replace MAX_KMER_SIZE in Makefile with an appropriate number or directly to make with

% make MAX_KMER_SIZE=64

In this case the maximum k-mer size allowed is 63, and each k-mer will use 16 bytes of memory.

Running
=======
To run the program use

% ./BFGraph

To list available commands

% ./BFGraph filter
% ./BFGraph contigs

will provide more complete help.

Documentation
=============
See docs/BFGDocumentation.pdf
You can recompile the documentation by running make inside the directory docs.

Example of usage
================
The directory example contains two small read files. 
See the file 'example.sh' to see how to run the program on these files.

Doxygen automatic documentation
===============================
Run `doxygen Doxyfile` to get html, latex and man documentation.

Notes
=====

* BFGraph was developed on x86-64 GNU/Linux. Porting to other unix-like
  platforms should be easy, but we haven't done so yet.

* If you run into bugs or problems or have suggestions for future versions 
  please contact me at pmelsted@gmail.com

License
=======


* The hash functions used are from the MurmurHash Library, version 3, released under the
  MIT License. http://code.google.com/p/smhasher/

* The kseq functions for reading fast(a|q)(.gz) files are copyrighted by Heng Li and released
  under the MIT License. http://lh3lh3.users.sourceforge.net/kseq.shtml


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
