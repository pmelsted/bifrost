# Bifrost API

Welcome to the introduction of the Bifrost API. In this tutorial, you will learn how to use the Bifrost C++ API to index, modify, query and color compacted de Bruijn graphs. Suggestions are always welcome so do not hesitate to leave some feedbacks.

### Before starting

This introduction assumes you are familiar with C++ although no advanced knowledge of the language is required. This tutorial is in construction so I will detail and clarify some parts of it over time. It is additionally assumed you are familiar with the general concept of de Bruijn graphs.

I will try to cover as much ground as possible in this tutorial. However, the keyword of this tutorial is 'simplicity' so I will not detail every single function of the API. If you are particularly interested in a function, the documentation (located in `/doc/doxygen/`) describes in detail what it does, what is the input and what is the output. If a function or a type requires special attention, I will point out that having a look at the documentation would be beneficial.

## Bifrost graphs

Graphs in Bifrost are compacted bi-directed de Bruijn graphs. It is important to understand this concept before getting started as the API relies on these notions. Unitigs, compaction and bi-directed edges are all explained in this excellent [short tutorial](https://github.com/GATB/bcalm/blob/master/bidirected-graphs-in-bcalm2/bidirected-graphs-in-bcalm2.md) by Paul Medvedev and Rayan Chikhi.
I will try to summarize quickly some of the important concepts for Bifrost:

- Bifrost uses a node-centric representation of the graph. Vertices are represented explicitely while edges are represented implicitly.
- Vertices are unitigs, i.e., sequences of length greater or equal to *k* (the *k*-mer size). Hence, a unitig is composed of at least one *k*-mer.
- *K*-mers in a unitig are not branching (edge in-degree = edge out-degree = 1) in the graph **except** the first and last *k*-mers which might be branching.
- A *k*-mer occurs in one unitig at most.
- Because the graph is bi-directed, a unitig represents two sequences: itself and its reverse-complement. Hence, a unitig can be traversed/read in two different ways: from left to right (*foward* or *+* direction) or from right to left by complementing the DNA symbols (*reverse* or *-* direction). This is the *strandness* of the unitig.
- Edges are directed: they have a start unitig *A* and an end unitig *B*. Furthermore, edges embed the strandness of the unitigs they connect, which can be: {+,+}, {+,-}, {-,+} and {-,-}. An edge *e* connecting the forward sequence of *A* to the reverse-complement sequence of *B* would be *e={A+,B-}*.

## Compacted de Bruijn graph

Let's start by creating a compacted de Bruijn graph.
```cpp
#include <bifrost/CompactedDBG.hpp>

int main(int argc, char **argv){

	CompactedDBG<> cdbg;
}
```

The first line includes the Bifrost C++ header for compacted de Bruijn graphs and that will be the only `#include` you need. Then, you just created an empty compacted de Bruijn graph named `cdbg` in the main function. That's it. Without parameters, that graph has a default *k*-mer size of `k=31` but you can change that later. Note that if you want to use a *k*-mer size larger than 31, you have to compile the Bifrost library with an extra parameter as explained [here](https://github.com/pmelsted/bifrost#large-k-mers).
The `<>` in `CompactedDBG<>` is a template for the type of data that you want to associate to unitigs. For now, we assume you do not want to associate data to unitigs so it is empty. Later in this tutorial, we will cover how to create a graph of type `CompactedDBG<MyData>` where type `MyData` represents some data you want to associate to unitigs.

In the next code samples, I will skip the include and the main function for simplicity, unless it is required to understand the code snippet.

Let's initialize our de Bruijn graph with a specific *k*-mer size that we will subsequently print.
```cpp
const size_t k = 31;

CompactedDBG<> cdbg(k);

cout << "K-mer size is " << cdbg.getK() << endl;
```

### Adding sequences

Now that we have an empty graph, let's add some data in it:
```cpp
const string seq = "ACGTCGTACGTCCCGTAAACGTTAAACGTAAACGTGTGTGCAAAATGTCTAGTTTTTTTACGCTGATATAGTC";

cdbg.add(seq);
```

Function `CompactedDBG::add()` takes as input a string (a DNA sequence) and inserts it in the graph. The function takes care of splitting the input sequence into unitigs. *K*-mers overlapping non-DNA characters will be automatically discarded.

Since it is the first time we modify a graph, here are two properties of Bifrost to keep in mind:
- Bifrost graphs are always compacted, no matter what. If you edit the graph by adding sequences, removing unitigs or merging graphs, Bifrost will *always* take care of compacting the graph. It is not possible to have an intermediate state where a Bifrost graph is not compacted.
- Most functions in the Bifrost API have an optional parameter *verbose* for printing information messages about the execution of the function. By default, `verbose=false` so no information messages are printed. By setting it to `true`, information message will be printed to the standard output `stdout`.

## Finding a ***k***-mer

Searching *k*-mers in the graph is one of the most common task to perform. It works as follows:
```cpp
const string kmer_sequence = "ACGTCGTACGTCCCGTAAACGTTAAACGTAA";
const Kmer km = Kmer(kmer_sequence.c_str());

UnitigMap<> um = cdbg.find(km);
```

An object `Kmer` is first created from the *k*-mer sequence and it is passed to function `CompactedDBG::find()` which returns a `UnitigMap` object. Objects of type `UnitigMap` are fundamental in Bifrost as they are used for nearly everything. A `UnitigMap` object represents the mapping of a sequence on a unitig of the graph. Hence, when searching for a *k*-mer in the graph, the returned `UnitigMap` object is the mapping of the searched *k*-mer sequence on a unitig of the graph or an empty mapping if the *k*-mer is not found:
```cpp
if (um.isEmpty) cout << "Kmer " << kmer_sequence << " was not found" << endl;
else {

	const string unitig = um.referenceUnitigToString();
	const string strandness = um.strand ? "forward" : "reverse-complement";
	const size_t position = um.pos;

	cout << "Kmer " << kmer_sequence << " was found in the " << strandness << " direction of unitig " << unitig << endl;
}
```

The important parameters of `UnitigMap` objects are:
- *isEmpty*: if `true`, the mapping is empty: the searched sequence was not found in the graph.
- *pos*: start position of the searched sequence on the unitig
- *len*: length of the mapping **in *k*-mers** (at least 1)
- *strand*: strandness of the mapped sequence, `true` for forward and `false` for reverse-complement
- *size*: length of the unitig **in bp** (at least *k*)

The type `UnitigMap` has tons of useful functions so it is probably a good idea to have a look at its documentation. Here is a quick overview:
- `UnitigMap::referenceUnitigToString()`: outputs the string of the unitig associated with this mapping.
- `UnitigMap::mappedSequenceToString()`: outputs the string of the mapped sequence. Takes into account the strandness.
- `UnitigMap::getUnitigHead()`: returns the first *k*-mer of the unitig associated with this mapping.
- `UnitigMap::getMappedHead()`: returns the first *k*-mer of the unitig in the mapped sequence. takes into account the mapping strandness.
- ...

As you can see above, some functions take into account the strandness and some functions do not. Here is a little bit of terminology to figure out what is what:
- reference unitig: the unitig to which the sequence map to, as stored in the graph (in forward direction)
- mapped sequence: a substring of the reference unitig, in forward or reverse-complement direction

All `UnitigMap` functions containing `referenceUnitig` or just `Unitig` return something related to the reference unitig. Hence, parameter `UnitigMap:strand` is not used here. All `UnitigMap` functions containing `mappedSequence` or just `Mapped` return something related to the mapped sequence (a substring of the reference unitig) and takes into account parameter `UnitigMap::strand`. Here is a figure to explain it all:

## Deleting sequences

Now that type `UnitigMap` has been introduced, we can keep going with how to delete a unitig from the graph. Let assume that we want to delete the unitig where was previously found our *k*-mer `km`. That unitig was the reference unitig of a `UnitigMap` object `um`. Removing that unitig works as follows:
```cpp
cdbg.remove(um);
```

It is as simple as that. Now, this function removes the unitig entirely. What if you want to remove only a substring, say just the *k*mer `km`?
```cpp
const string unitig = um.referenceUnitigToString();
const string unitig_pre = unitig.substr(0, um.pos + k - 1); // Unitig substring 'before' k-mer
const string unitig_suf = unitig.substr(um.pos + 1, um.size - um.pos - 1); // Unitig substring 'after' k-mer

cdbg.remove(um);
cdbg.add(unitig_pre);
cdbg.add(unitig_suf);
```

To remove a substring, you must remove the unitig entirely and re-insert the unitig parts that came before and after the substring.

## Building the graph from a dataset

## Storing unitigs identity

A question that often comes back to me is 'How do I save a unitig identity or its position in the graph?'. In C++, if you have a vector of objects and you want to remember a specific object for later, you store an iterator to that object or the position of that object in the vector. Let assume that for the purpose of your program, you want to save the reference unitig of a `UnitigMap` object `um` where you found your *k*-mer `km`:
```cpp
vector<UnitigMap<>> v_um;

v_um.push_back(um);
```

That is the closest equivalent of storing an iterator, quick and fast. However, beware that same as for a vector iterator, if you decide to modify the graph (adding or removing a sequence) after that, the stored `UnitigMap` will be broken and it will not point anymore to a valid location in the graph. Using the stored `UnitigMap` object after modifying the graph is undefined behavior and your program might very well crash. The reason for that is that when you modify the graph, it is automatically re-compacted: some unitigs might split or merge. Hence, the unitig pointed out by the stored `UnitigMap` object might look different or might not exist at all anymore.

There is another way to store unitig identity or positions which will be ideal in a number of situations. Remember that each *k*-mer in the graph occur in **at most** one unitig. Which means that a *k*-mer can be used as an identifier for a unitig:
```cpp
vector<Kmer> v_km;

for (auto& um : v_um) {

	const Kmer head_kmer = um.getUnitigHead(); // First k-mer of the unitig

	v_km.push_back(head_kmer); // Push the head k-mer to a vector
}

v_um.clear(); // Vector of UnitigMap is not needed anymore
```
Here, we use the head *k*-mer of the unitig as its identifier. Hence, instead of having a vector of `UnitigMap`, you have a vector of `Kmer` where each *k*-mer is the head *k*-mer of a unitig you are interested in. Now to retrieve the unitigs associated to these *k*mers:
```cpp
for (auto& km : v_km) {

	UnitigMap<> um = cdbg.find(km, true);
}
```

The major advantage of using `Kmer` rather than `UnitigMap` is that a `Kmer` object is a lot less memory consuming than storing a `UnitigMap` object, which is useful for long vectors of unitigs. On the other hand, retrieving the unitig associated to a `Kmer` costs some time. It is a tradeoff you have to decide for yourself. My advise: if you have only a few hundreds/thousands of unitigs to remember, store `UnitigMap` objects. Otherwise, store `Kmer` objects.

Not that in the previous code snippet, I used `cdbg.find(km, true)`. When this last parameter is set to `true (`false` by default), it indicates that you want to search for this *k*-mer **only** at the extremities of unitigs (head or tail *k*-mers only). Doing this significantly speeds up the search and it comes very handy when you search for *k*-mers that you know are the head *k*-mer of unitigs.