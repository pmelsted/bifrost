#include <cstdio>
#include <ctime>
#include <iostream>
#include <map>
#include <cassert>

#include "../CompressedCoverage.hpp"

using namespace std;
// strings created in python with:
// f = lambda x,start,new: x[0:start*2] + new + x[start*2+len(new):]
// x = "00"*28
// f(x,12,"01"*12)
//
int main(int argc, char *argv[]) {
  CompressedCoverage c1(28);
  c1.cover(12,24);
  string s1 = "00000000000000000000000001010101010101010101010100000000";
  string s1r = c1.toString();
  printf("toString: %s\n", c1.toString().c_str());
  cout << "s1: " << s1 << endl;
  for (size_t index = 0; index < 30; ++index) {
    printf("sir[%u] = %c\n", 2*index+1, s1r[2*index+1]);
    //assert(s1r[2*index+1] == '1');
  }
  assert(c1.toString() == s1);
  string s2 = "0000000000000000000000000101010101010101010101010000000000";
  CompressedCoverage c2(29);
  c2.cover(12,24);
  assert(c2.toString() == s2);

  cout << &argv[0][2] << " completed successfully" << endl;
}
