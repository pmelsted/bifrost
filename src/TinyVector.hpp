#ifndef KALLISTO_TINYVECTOR_H
#define KALLISTO_TINYVECTOR_H

#include <vector>
#include <stdint.h>

template<class T, int N = (1 + (48 / sizeof(T)))>
class tiny_vector {
public:
  tiny_vector() : short_(true) {arr.size = 0;};
  tiny_vector(std::initializer_list<T> l) {
    short_ = true;
    arr.size = 0;
    reserve(l.size());
    for (auto &x : l) {
      push_back(x);
    }
  }

  ~tiny_vector() {
    clear();
  };
  union {
    struct {
      T data[N];
      uint8_t size;
    } arr;
    struct {
      T *data;
      size_t cap;
      size_t size; // 24 bytes + ~24 malloc overhead
    } vec;
  };
  tiny_vector(const tiny_vector& o) {
    short_ = true;
    arr.size = 0;
    reserve(o.capacity());
    for (auto &x : o) {
      push_back(x);
    }
  }
  tiny_vector(const std::vector<T>& o) {
    short_ = true;
    arr.size = 0;
    reserve(o.capacity());
    for (auto &x : o) {
      push_back(x);
    }
  }
  tiny_vector(tiny_vector&& o) {
    if (o.isShort()) {
      size_t sz = o.size();
      for (size_t i = 0; i < sz; i++) {
        arr.data[i] = std::move(o.arr.data[i]);
      }
      o.arr.size = 0;
      arr.size = sz;
      short_ = true;
    } else {
      vec.data = o.vec.data;
      vec.size = o.vec.size;
      vec.cap  = o.vec.cap;
      o.vec.data = nullptr;
      o.arr.size = 0;
      o.short_ = true;
      short_ = false;
    }
  }
  tiny_vector& operator=(const tiny_vector& o) {
    clear();
    reserve(o.capacity());
    for (auto &x : o) {
      push_back(x);
    }
    return *this;
  }

  bool short_;

  typedef T* iterator;
  typedef T const* const_iterator;


  iterator begin() { return getPointer(); }
  const_iterator begin() const { return getPointer(); }
  iterator end() { return getPointer() + size(); }
  const_iterator end() const { return (getPointer() + size()); }


  bool isShort() const {
    return short_;
  }



  T* getPointer() {
    if (!isShort()) {
      return vec.data;
    } else {
      return &arr.data[0];
    }
  }

  const T* getPointer() const {
    if (!isShort()) {
      return vec.data;
    } else {
      return &arr.data[0];
    }
  }

  size_t size() const {
    if (isShort()) {
      return arr.size;
    } else {
      return vec.size;
    }
  };

  bool empty() const {
    return (size()==0);
  }

  bool operator==(const tiny_vector& o) const {
    if (size() == o.size()) {
      for (auto it=begin(), oit = o.begin(); it != end(); ++it,++oit) {
        if (*it == *oit) {
          // pass
        } else {
          return false;
        }
      }
      return true;
    } else {
      return false;
    }
  }

  bool operator!=(const tiny_vector& o) const {
    return !operator==(o);
  }

  size_t capacity() const {
    if (isShort()) {
      return N;
    } else {
      return vec.cap;
    }
  }

  T& operator[](size_t i) {
    return *(begin() + i);
  };

  const T& operator[](size_t i) const {
    return *(begin() + i);
  }

  void push_back(const T& value) {
    if (size() < capacity()) {
      // we have room for one more value
    } else {
      _reserve_and_copy(3*size()/2);
    }
    *(end()) = value;
    //new(end()) T(std::forward(value)); // fancy way of assigning to pointer
    if (isShort()) {
      ++arr.size;
    } else {
      ++vec.size;
    }
  };

  //void push_back(T&& value) {}


  void clear() {
    if (!isShort()) {
      delete[] vec.data;
      vec.data = nullptr;
      vec.size = 0;
      vec.cap = 0;
    } else {
      for (auto it = begin(); it != end(); ++it) {
        it->~T(); // manually call destructor
      }
      arr.size = 0;
    }
    short_ = true;
  }

  void reserve(size_t sz) {
    if (sz > capacity()) {
      _reserve_and_copy(sz);
    }
  }

  void _reserve_and_copy(size_t sz) {
    if (sz <= capacity() ) {
      return;
    }
    size_t old_size = size();

    T* newdata = new T[sz];
    // move old elements over
    size_t i = 0;
    for (auto it = begin(); it != end(); ++it,++i) {
      newdata[i] = std::move(*it); //
    }
    if (!isShort()) {
      delete[] vec.data;
    }

    short_ = false;
    vec.size = old_size;
    vec.cap = sz;
    vec.data = newdata;
  }


};

#endif
