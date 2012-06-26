#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <omp.h>

using namespace std;

int main () {
  int MAX = 100, reads = 0;
  int l[MAX];
  int sum = 0;
  vector<int> *p;
  vector<int> *v;
  size_t threadcount = 0;

  while (reads < MAX && cin >> l[reads]) {
    ++reads;
  }

  #pragma omp parallel
  {
    #pragma omp master 
    {
      threadcount = omp_get_num_threads();
      cout << "Number of threads : " << threadcount << endl;
    }
  }
  v = new vector<int>[threadcount];

  #pragma omp parallel default(shared) private(p)
  {
    p = &v[omp_get_thread_num()];
    #pragma omp for nowait
    for (int i=0; i<reads; i++) {
      printf("l[%d]=%d\n", i, l[i]);
      sum += l[i];
      p->push_back(l[i]);
    }
  }  


  #pragma omp parallel shared(sum)
  {
    cout << "I'm a thread!\n";
  }

  cout << "Sum: " << sum << endl;
  for (int i=0; i<threadcount; i++) {
    for(vector<int>::iterator it=v[i].begin(); it != v[i].end(); ++it) {
      printf("from vector: %d\n", *it);
    }
  }

  delete[] v;
}
