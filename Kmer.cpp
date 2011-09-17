
#include "Kmer.hpp"

char *int2bin(uint32_t a, char *buffer, int buf_size) {
  //buffer += (buf_size - 1);

    for (int i = 7; i >= 0; i--) {
        *buffer++ = (a & 1) + '0';
        a >>= 1;
    }

    return buffer;
}

static const uint8_t base_swap[256] = {
	0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0,
	0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70, 0xb0, 0xf0,
	0x04, 0x44, 0x84, 0xc4, 0x14, 0x54, 0x94, 0xd4,
	0x24, 0x64, 0xa4, 0xe4, 0x34, 0x74, 0xb4, 0xf4,
	0x08, 0x48, 0x88, 0xc8, 0x18, 0x58, 0x98, 0xd8,
	0x28, 0x68, 0xa8, 0xe8, 0x38, 0x78, 0xb8, 0xf8,
	0x0c, 0x4c, 0x8c, 0xcc, 0x1c, 0x5c, 0x9c, 0xdc,
	0x2c, 0x6c, 0xac, 0xec, 0x3c, 0x7c, 0xbc, 0xfc,
	0x01, 0x41, 0x81, 0xc1, 0x11, 0x51, 0x91, 0xd1,
	0x21, 0x61, 0xa1, 0xe1, 0x31, 0x71, 0xb1, 0xf1,
	0x05, 0x45, 0x85, 0xc5, 0x15, 0x55, 0x95, 0xd5,
	0x25, 0x65, 0xa5, 0xe5, 0x35, 0x75, 0xb5, 0xf5,
	0x09, 0x49, 0x89, 0xc9, 0x19, 0x59, 0x99, 0xd9,
	0x29, 0x69, 0xa9, 0xe9, 0x39, 0x79, 0xb9, 0xf9,
	0x0d, 0x4d, 0x8d, 0xcd, 0x1d, 0x5d, 0x9d, 0xdd,
	0x2d, 0x6d, 0xad, 0xed, 0x3d, 0x7d, 0xbd, 0xfd,
	0x02, 0x42, 0x82, 0xc2, 0x12, 0x52, 0x92, 0xd2,
	0x22, 0x62, 0xa2, 0xe2, 0x32, 0x72, 0xb2, 0xf2,
	0x06, 0x46, 0x86, 0xc6, 0x16, 0x56, 0x96, 0xd6,
	0x26, 0x66, 0xa6, 0xe6, 0x36, 0x76, 0xb6, 0xf6,
	0x0a, 0x4a, 0x8a, 0xca, 0x1a, 0x5a, 0x9a, 0xda,
	0x2a, 0x6a, 0xaa, 0xea, 0x3a, 0x7a, 0xba, 0xfa,
	0x0e, 0x4e, 0x8e, 0xce, 0x1e, 0x5e, 0x9e, 0xde,
	0x2e, 0x6e, 0xae, 0xee, 0x3e, 0x7e, 0xbe, 0xfe,
	0x03, 0x43, 0x83, 0xc3, 0x13, 0x53, 0x93, 0xd3,
	0x23, 0x63, 0xa3, 0xe3, 0x33, 0x73, 0xb3, 0xf3,
	0x07, 0x47, 0x87, 0xc7, 0x17, 0x57, 0x97, 0xd7,
	0x27, 0x67, 0xa7, 0xe7, 0x37, 0x77, 0xb7, 0xf7,
	0x0b, 0x4b, 0x8b, 0xcb, 0x1b, 0x5b, 0x9b, 0xdb,
	0x2b, 0x6b, 0xab, 0xeb, 0x3b, 0x7b, 0xbb, 0xfb,
	0x0f, 0x4f, 0x8f, 0xcf, 0x1f, 0x5f, 0x9f, 0xdf,
	0x2f, 0x6f, 0xaf, 0xef, 0x3f, 0x7f, 0xbf, 0xff
};

Kmer::Kmer() {
  memset(bytes,0,MAX_K/4);
}

Kmer::Kmer(const Kmer& o) {
  memcpy(bytes,o.bytes,MAX_K/4);
}

Kmer::Kmer(const char *s) {
  set_kmer(s);
}

Kmer& Kmer::operator=(const Kmer& o) {
  if (this != &o) {
    memcpy(bytes, o.bytes, MAX_K/4);
  }
  return *this;
}

void Kmer::set_deleted() {
  memset(bytes,0xff,MAX_K/4);
}

bool Kmer::operator<(const Kmer& o) const {

  return memcmp(bytes,o.bytes,MAX_K/4) < 0;
    /*
    //printf("Less\n");
        
    for (int i = 0; i < k_bytes-1; ++i) {
      //printf("  %x vs %x\n",bytes[i],o.bytes[i]);
      if (bytes[i] >= o.bytes[i]) {
        return false;
      }
    }

    //printf("first match, but we have some bytes");
    //printf(" %d vs %d\n",bytes[k_bytes-1],o.bytes[k_bytes-1]);
    //printf(" now with masks ");
    //printf(" %d vs %d\n",bytes[k_bytes-1] &k_modmask,o.bytes[k_bytes-1]&k_modmask);
    if ((bytes[k_bytes-1] & k_modmask) >= (o.bytes[k_bytes-1] & k_modmask)) {
      return false;
    }
    return true;*/
}

bool Kmer::operator==(const Kmer& o) const {
  
  return memcmp(bytes,o.bytes,MAX_K/4)==0;
  /*
  for (int i = 0; i < k_bytes-1; ++i) {
      if (bytes[i] != o.bytes[i]) {
        return false;
      }
    }
    if ((bytes[k_bytes-1] & k_modmask) != (o.bytes[k_bytes-1] & k_modmask)) {
      return false;
    }
    return true;
  */
}

void Kmer::set_kmer(const char *s)  {
  size_t i,j,l;
  //printf("%s\n",s);
  memset(bytes,0,MAX_K/4);

  //printf("before   "); this->printBinary();
    
  for (i = 0; i < k; ++i) {
    j = i % 4;
    l = i/4;
    assert(*s != '\0');
    switch(*s) {
    case 'A': break;
    case 'C': bytes[l] |= (0x01 << (2*j)); break;
    case 'G': bytes[l] |= (0x02 << (2*j)); break;
    case 'T': bytes[l] |= (0x03 << (2*j)); break;
    }
    //printf("step %2d %c",i,*s);  printBinary();
      
    s++; 
  }
}

uint64_t Kmer::hash() const {
  //this->printBinary();
  //for (int i = 0; i < k_bytes; ++i) {
  //printf("%02X ",bytes[i]);
  //}
  //printf("\n");
  uint64_t ret;
  MurmurHash3_x64_64((const void*)bytes,k_bytes,0,&ret);
  return ret;
}


Kmer Kmer::twin() const {
    
  Kmer km(*this);
  // char s[1024];
    
//      int2bin(k_modmask,s,32);
//     printf("modmask: %s\n",s);
//     printf("TWIN:\n"); 
//     printf("original"); this->printBinary();
//     printf("before  "); km.printBinary();
    
  for (size_t i = 0; i < k_bytes; i++) {
    km.bytes[i] = ~bytes[i];
  }
  // printf("~       "); km.printBinary();
  km.bytes[k_bytes-1] ^= ~k_modmask;
  // printf("!       ");     km.printBinary();
    
  // shift to byte alignment
  km.shiftLeft(8*k_bytes-2*k);

  //printf("shift   ");    km.printBinary();

  uint8_t tmp;
  for (size_t i = 0; i < k_bytes/2; ++i) {
    tmp = km.bytes[i];
    km.bytes[i] = base_swap[km.bytes[k_bytes-1-i]];
    km.bytes[k_bytes-1-i] = base_swap[tmp];
  }

  if ((k_bytes %2) == 1) {
    km.bytes[k_bytes/2] = base_swap[km.bytes[k_bytes/2]];
  }

  //  printf("rotated ");    km.printBinary(); printf("\n");

  return km;
}


Kmer Kmer::getLink(const size_t index) const {
  assert(index >=0 && index < 8);
  size_t ind_mod = index % 4;
  char c = '\0';
  switch (ind_mod) {
  case 0: c = 'A'; break;
  case 1: c = 'C'; break;
  case 2: c = 'G'; break;
  case 3: c = 'T'; break;
  }
  
  if (index < 4) {
    return forwardBase(c);
  } else {
    return backwardBase(c);
  }
}



Kmer Kmer::forwardBase(const char b) const {
  //  char tmp[1024];
  //this->toString(tmp);
  //printf("\n\nforward\n%s\n",tmp);
  // printf("adding %c\n",b);

  int s = 2*((k+3) % 4);
  //printf("shift used %d\n",s);

  Kmer km(*this);
  //printf("before "); km.printBinary();

  km.shiftRight(2);
  //printf("shift  "); km.printBinary();
  km.bytes[k_bytes-1] &= Kmer::k_modmask;

  //printf("mask   "); km.printBinary();
  switch(b) {
  case 'A': km.bytes[k_bytes-1] |= 0x00 << s; break;
  case 'C': km.bytes[k_bytes-1] |= 0x01 << s; break;
  case 'G': km.bytes[k_bytes-1] |= 0x02 << s; break;
  case 'T': km.bytes[k_bytes-1] |= 0x03 << s; break;
  }
  //printf("after  "); km.printBinary();
  //km.toString(tmp);
  //printf("%s\n\n",tmp);
  return km;
}

Kmer Kmer::backwardBase(const char b) const {
  // char tmp[1024];
//   this->toString(tmp);
//   printf("backward\n%s\n",tmp);
//   printf("adding %c\n",b);

  //size_t s= 2*((k+3)%4);
  //printf("shift used %d\n",s);

  Kmer km(*this);
  //printf("before "); km.printBinary();

  km.shiftLeft(2);
  // printf("shift  "); km.printBinary();
  km.bytes[k_bytes-1] &= Kmer::k_modmask;
  if (k%4 == 0 and k_bytes < MAX_K/4) {
    km.bytes[k_bytes] = 0x00;
  }
  //printf("mask   "); km.printBinary();

  switch(b) {
  case 'A': km.bytes[0] |= 0x00; break;
  case 'C': km.bytes[0] |= 0x01; break;
  case 'G': km.bytes[0] |= 0x02; break;
  case 'T': km.bytes[0] |= 0x03; break;
  }
  // printf("after  "); km.printBinary();
  // km.toString(tmp);
  // printf("%s\n",tmp);

  return km;
}

void Kmer::printBinary() const {
  char buff[9]; buff[8] = '\0';
  printf("binary:");
  for (size_t i = 0; i < Kmer::k_bytes; i++) {
    int2bin(bytes[i],buff,8);
    printf("%s",buff);
  }
  printf("\n");
}

void Kmer::toString(char * s) const {
  size_t i,j,l;
  
  for (i = 0; i < k; i++) {
    j = i % 4;
    l = i / 4;
    switch(((bytes[l]) >> (2*j) )& 0x03 ) {
    case 0x00: *s = 'A'; ++s; break;
    case 0x01: *s = 'C'; ++s; break;
    case 0x02: *s = 'G'; ++s; break;
    case 0x03: *s = 'T'; ++s; break;  
    }
  }
  *s = '\0';

    /*
      printf("bytes: ");
    for (i = 0; i < k_bytes; i++ ) {
      printf("%x ",bytes[i]);
    }
    printf("\n");
    */

}


void Kmer::shiftLeft(int shift) {
  if (shift>0) {
    if (shift < 8 ) {
      for (size_t i = Kmer::k_bytes-1; i > 0; i--) {
	bytes[i] <<= shift;
	bytes[i] |= (uint8_t) (bytes[i-1] >> (8-shift));
      }
      bytes[0] <<= shift;
    } else {
      // we should never need this!
      assert(0);
    }
  }
}

void Kmer::shiftRight(int shift) {
  if (shift >0) {
    if (shift < 8) {
      for (size_t i = 0; i < Kmer::k_bytes-1; i++) {
	bytes[i] >>= shift;
	bytes[i] |= (uint8_t) ( bytes[i+1] << (8-shift));
      }
      bytes[Kmer::k_bytes-1] >>= shift;
    } else {
      // bad
      assert(0);
    }
  }
}

void Kmer::set_k(unsigned int _k) {
  assert(_k <= MAX_K);
  assert(_k > 0);
  assert(k_bytes == 0); // we can only call this once
  k = _k;
  k_bytes = (_k+3)/4;
  //  k_longs = (_k+15)/16;
  k_modmask = (1 << (2*((k%4)?k%4:4)) )-1;
  
  //printf("k: %d, k_bytes: %d\n",k,k_bytes);
}



unsigned int Kmer::k = 0;
unsigned int Kmer::k_bytes = 0;
//unsigned int Kmer::k_longs = 0;
unsigned int Kmer::k_modmask = 0;
