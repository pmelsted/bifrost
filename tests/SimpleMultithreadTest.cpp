#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>

#include "../MutexLock.hpp"

using namespace std;

int main () {
  int MAX = 100, counter = 0, reads = 0;
  int l[MAX];
  int index = 0;
  int sum = 0;
  vector<int> v, vi;

  while (reads < MAX && cin >> l[index++]) {
    ++reads;
  }
  
  MutexLock ml;

  #pragma omp parallel default(shared) private(vi)
  {
    #pragma omp for nowait
    for (int i=0; i<reads; i++) {
      sum += l[i];
      vi.push_back(l[i]);
    }

    ml.lock(); 
    v.insert(v.end(), vi.begin(), vi.end()); 
    ml.unlock();
  }  


  #pragma omp parallel shared(counter)
  {
    ++counter;
  }
  cout << "Number of cores: " << counter << endl;
  cout << "Sum: " << sum << endl;

  for(vector<int>::iterator it=v.begin(); it != v.end(); ++it) {
    printf("from vector: %d\n", *it);
  }
}
