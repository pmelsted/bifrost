#ifdef _OPENMP
#include <omp.h>
struct MutexLock
{
  MutexLock() { omp_init_lock(&innerlock); }
  ~MutexLock() { omp_destroy_lock(&innerlock); }
  void lock() { omp_set_lock(&innerlock); }
  void unlock() { omp_unset_lock(&innerlock); }
 
  MutexLock(const MutexLock& ) { omp_init_lock(&innerlock); }
  MutexLock& operator= (const MutexLock& ) { return *this; }
public:
  omp_lock_t innerlock;
};
#else
struct MutexLock
{
  void Lock() {}
  void Unlock() {}
};
#endif
