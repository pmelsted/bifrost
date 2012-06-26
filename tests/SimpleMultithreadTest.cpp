#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>

using namespace std;

int main () {
  int MAX = 100, reads = 0;
  int l[MAX];
  int sum = 0;
  vector<int> v, vi;

  while (reads < MAX && cin >> l[reads++]);

  #pragma omp parallel default(shared) private(vi)
  {
    #pragma omp for nowait
    for (int i=0; i<reads; i++) {
      sum += l[i];
      vi.push_back(l[i]);
    }

    #pragma omp master 
    {
      cout << "Number of threads : " << omp_get_num_threads() << endl;
    }

    #pragma omp critical 
    { 
      v.insert(v.end(), vi.begin(), vi.end()); 
    }
  }  


  #pragma omp parallel shared(sum)
  {
    cout << "I'm a thread!" << endl;
  }

  cout << "Sum: " << sum << endl;

  for(vector<int>::iterator it=v.begin(); it != v.end(); ++it) {
    printf("from vector: %d\n", *it);
  }
}
