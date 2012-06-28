#ifndef BFG_COMPRESSED_COVERAGE_HPP
#define BFG_COMPRESSED_COVERAGE_HPP

#include <stdint.h>
#include <cstring>


/* Short description: 
 *  - Tagged pointer union that is either
 *    - a pointer to a char array that stores 2-bit integers
 *    - a local 2-bit integer array
 *  - The bits are stored either as
 *    pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|pppppppp|ppppppF0
 *    dddddddd|dddddddd|dddddddd|dddddddd|dddddddd|dddddddd|dddddddd|ssssssF1
 *
 *    - First bit is 0 for pointer and 1 for local array
 *    - Second bit is 1 for a full coverage and 0 for non-full
 *    - For local array next 6 bits store size of the array
 *    - For the pointer version next 62 bits encode the pointer (last two bits
 *      are zeroed out before dereferencing)
 *    - The pointer points to an array of bytes, where the first 8 encode the
 *      size of the array used in uint32_t and the number of full positions.
 *    - The remainder of bytes are 2-bit encoded integers.
 *    - If the full bit is set then the pointer must be 0 and the memory released
 * */
class CompressedCoverage {
public: 
  CompressedCoverage(size_t size);
  ~CompressedCoverage();

  void cover(size_t start, size_t end); 

  size_t size() const;
  std::string toString() const; // for debugging

private:

  static const intptr_t tagMask = 1; // local array bit
  static const intptr_t fullMask = 2; // full bit
  static const intptr_t sizeMask = 0xFC; // 0b11111100
  static const intptr_t pointerMask = ~(tagMask | fullMask); // rest of bits

  static const size_t size_limit = 28; // 56 bit array, 28 2-bit integers
  
  uint8_t* getPointer() const;
  bool isFull() const;
  void releasePointer();
  
  union {
    uint8_t *asPointer;
    intptr_t asBits;
  };
};

#endif // BFG_COMPRESSED_COVERAGE_H
