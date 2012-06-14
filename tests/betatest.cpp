#include <cassert>

static const char alpha[4] = {'A', 'C', 'G', 'T'};
static const char beta[5] = {'T','G','A','C','N'}; 

char f(char c) {
  return ((c << 3) % 5) & 7;
}

char g(char c) {
  return (x & 7) >> 1;
}

int main(void) {
  assert(beta[f('A')] == 'T');
  assert(beta[f('C')] == 'G');
  assert(beta[f('G')] == 'C');
  assert(beta[f('T')] == 'A');
  assert(beta[f('N')] == 'N');
  
  assert(beta[g('A')] == 'T');
  assert(beta[g('C')] == 'G');
  assert(beta[g('G')] == 'C');
  assert(beta[g('T')] == 'A');
}

