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

// Not done:
// 2* 1 * 2
// 2* 1 * 2
// 2* 1* 2* 1* 2*
int k;

int test_checkEndKmer() {
  return 0;
}

// test checkEndKmer
// why not work  ? some infiniteloop
int test_splitAllContigs() {
  string s("CGCCGCCATGCCGACCATCCCTTTCATCCCCGTACCAGACACGCTGACCATTGCCATGTTATTCAGATTGCATTTATTGAGAAATTCTCCTTTTTCCCTG");
  cout << "Covering: " << s << endl;
  assert(s.size() == 100);
  BloomFilter bf;
  ContigMapper cmap;
  cmap.mapBloomFilter(&bf);
  Contig *c = new Contig(s.c_str());
  Kmer head(s.c_str());

  c->ccov.cover(0,69); // covers kmers [0,69]
  for (int i = 0; i <= 69; i++) {
    assert(c->ccov.covAt(i) == 1);
  }

  c->ccov.cover(0,69); // covers [0,69]

  /* We cannot do this because covAt asserts that the coverage isn't full
  for (int i = 0; i <= 69; i++) {
    assert(c->ccov.covAt(i) == 2);
  }
  */

  cmap.lContigs.insert(make_pair(head, c));
  const char* cs = s.c_str();
  size_t len = s.size();

  // we must insert the shortcuts as well!
  for (size_t i = k; i < len-k; i += k) {
    cmap.shortcuts.insert(make_pair(Kmer(cs+i),make_pair(head,i)));
  }
  // also insert shortcut for last k-mer
  cmap.shortcuts.insert(make_pair(Kmer(cs+len-k),make_pair(head,len-k)));

  pair<size_t, size_t> p = cmap.splitAllContigs(); // splited, deleted
  cout << p.first << " " << p.second << endl;

  // We have full coverage so nothing should be splitted, and nothing deleted
  assert(p.first == 0); // splitted
  assert(p.second == 0); // deleted
  //      delete c; // who deletes this contig? // double free :/ 
  return 0;
}

// 1 1* 2*1* 1
int test_splitAllContigs121() {
  string s("CGCCGCCATGCCGACCATCCCTTTCATCCCCGTACCAGACACGCTGACCATTGCCATGTTATTCAGATTGCATTTATTGAGAAATTCTCCTTTTTCCCTG");
  assert(s.size() == 100);
  // Kmer0 [0,31[
  // Kmer1 [1,32[
  // ...
  // Kmer9 [10,41[
  // Kmer69 [70,101[
  cout << s.size() << endl;
  for (int i=0; i <= 69; i++) {
    for (int j=i+1; j < 69; j++) {
      printf("(i,j)=(%d,%d)\n", i,j);
      BloomFilter bf;
      ContigMapper cmap;
      cmap.mapBloomFilter(&bf);
      Contig *c = new Contig(s.c_str());
      Kmer head(s.c_str());
      c->ccov.cover(0,69); // covers [0,69]
      c->ccov.cover(i,j);
      cmap.lContigs.insert(make_pair(head, c));
      const char* cs = s.c_str();
      size_t len = s.size();

      // we must insert the shortcuts as well!
      for (size_t i = k; i < len-k; i += k) {
        cmap.shortcuts.insert(make_pair(Kmer(cs+i),make_pair(head,i)));
      }
      // also insert shortcut for last k-mer
      cmap.shortcuts.insert(make_pair(Kmer(cs+len-k),make_pair(head,len-k)));
      pair<size_t, size_t> p = cmap.splitAllContigs();
      cout << p.first << " " << p.second << endl;
      assert(p.first == 1); // splitted
      assert(p.second == 0); // deleted
      //      delete c; // who deletes this contig? // double free :/ 
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
  k = 31;
  Kmer::set_k(k);
  test_splitAllContigs();
  test_splitAllContigs121();
  return 0;
}
