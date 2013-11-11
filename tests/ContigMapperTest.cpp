#include <cassert>

#define private public
#include "../ContigMapper.hpp"

// splitAllContigs 
//
// Búa til Contiga alveg sama hvernig tengdir.
// Með coverage 1 einhvers staðar
// Með coverage 2 alls stadar 
// Með coverage 1 alls stadar 
// Með coverage 0 alls stadar 
// asserta ad allt splittast rett
// checka ad head og tail mappist rett i nyja contignum
// Bua til cmap
// Bua til stutta og langa contiga
//

int test_checkEndKmer() {
  return 0;
}

// test checkEndKmer
// why not work  ???
int test_splitAllContigs() {
// 2* 1 * 2
// 2* 1* 2* 1* 2*
  //string s = string(100,'A');
  string s("CGCCGCCATGCCGACCATCCCTTTCATCCCCGTACCAGACACGCTGACCATTGCCATGTTATTCAGATTGCATTTATTGAGAAATTCTCCTTTTTCCCTG");
  assert(s.size() == 100);
  cout << s.size() << endl;
  for (int i=0; i < 100; i++) {
    for (int j=i+1; j < 100; j++) {
      // 1 2* 1
      // 2 on [i,j[
      BloomFilter bf;
      ContigMapper cmap;
      cmap.mapBloomFilter(&bf);
      Contig *c = new Contig(s.c_str());
      Kmer head(s.c_str());
      c->ccov.cover(0,100); // covers [0,100[
      c->ccov.cover(i,j);
cmap.lContigs.insert(make_pair(head, c));
//      delete c; // who deletes this contig? // double free :/ 
      pair<size_t, size_t> p = cmap.splitAllContigs();
      cout << p.first << " " << p.second << endl;
    }
  }

  
  cout << s << endl;
  return 0;
}

int test_fwBfStep() {
//  ContigMapper cmap;
//  BloomFilter bf
  return 0;

}

int main(int argc, char **argv) {
  test_splitAllContigs();
  return 0;
}
