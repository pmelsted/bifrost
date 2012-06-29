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
  string s1r = c1.toString();
  assert(s1r.substr(0, 28*2) == string(2*28, '0'));

  c1.cover(27,27);
  c1.cover(27,27);
  s1r = c1.toString();
  assert(s1r[54 - 27*2] == '1');      // 0
  assert(s1r[54 - 27*2 +1] == '0');   // 1
  
  c1.cover(0,0);
  s1r = c1.toString();
  assert(s1r[54 - 0*2] == '0');       // 54
  assert(s1r[54 - 0*2 +1] == '1');    // 55
  
  c1.cover(1,1);
  s1r = c1.toString();
  assert(s1r[54 - 1*2] == '0');      // 52
  assert(s1r[54 - 1*2 +1] == '1');   // 53

  c1.cover(10,14);
  c1.cover(10,14);
  c1.cover(10,14);
  c1.cover(10,14);
  s1r = c1.toString();
  for(size_t i = (54 - 10*2); i >= (54 - 14*2); i-=2) {
    assert(s1r[i] == '1');
    assert(s1r[i+1] == '0');
  }

  
  c1.cover(20,23);
  s1r = c1.toString();
  for(size_t i = (54 - 20*2); i >= (54 - 23*2); i-=2) {
    assert(s1r[i] == '0');
    assert(s1r[i+1] == '1');
  }
  
  CompressedCoverage c2(10);
  string s2r = c2.toString();
  assert(s2r.substr(0, 28*2) == string(2*28, '0'));

  c2.cover(9,9);
  c2.cover(9,9);
  s2r = c2.toString();
  assert(s2r[54 - 9*2] == '1');
  assert(s2r[54 - 9*2 +1] == '0');
  
  c2.cover(0,0);
  s2r = c2.toString();
  assert(s2r[54 - 0*2] == '0');
  assert(s2r[54 - 0*2 +1] == '1');
  
  c2.cover(1,1);
  s2r = c2.toString();
  assert(s2r[54 - 1*2] == '0');
  assert(s2r[54 - 1*2 +1] == '1');

  c2.cover(5,7);
  c2.cover(5,7);
  s2r = c2.toString();
  for(size_t i = (54 - 5*2); i >= (54 - 7*2); i-=2) {
    assert(s2r[i] == '1');
    assert(s2r[i+1] == '0');
  }


  CompressedCoverage c3(50);
  string s3r = c3.toString();
  cout << s3r << endl;
  c3.cover(10,20);
  c3.cover(30,40);
  s3r = c3.toString();
  cout << s3r << endl;

  cout << &argv[0][2] << " completed successfully" << endl;
}
