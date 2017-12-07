#include "BlockedBloomFilter.hpp"

#if defined(__AVX2__)

const __m256i BlockedBloomFilter::mask_and_div = _mm256_set1_epi16(MASK_BITS_BLOCK);
// All 4 LSB of each 16 bits word
const __m256i BlockedBloomFilter::mask_and_mod = _mm256_set1_epi16(0xf);
// All 1 LSB of each 32 bits word
const __m256i BlockedBloomFilter::one2shift_lsb = _mm256_set1_epi32(0x1);
// Set the 17th LSB bit of each 32 bits word
const __m256i BlockedBloomFilter::one2shift_msb = _mm256_set1_epi32(0x10000);
// All 16 LSB bits of each 32 bits word
const __m256i BlockedBloomFilter::mask_lsb = _mm256_set1_epi64x(0x0000ffff0000ffff);

#endif
