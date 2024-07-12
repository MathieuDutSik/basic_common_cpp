// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_HASH_FUNCTIONS_H_
#define SRC_BASIC_HASH_FUNCTIONS_H_

// clang-format off
#include "Temp_common.h"
#include <array>
#include <cstdint>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

// Code copy pasted from
// https://codereview.stackexchange.com/questions/163872/md5-implementation-in-c11

class md5 {
private:
  std::uint32_t a0_;
  std::uint32_t b0_;
  std::uint32_t c0_;
  std::uint32_t d0_;

  std::array<std::uint32_t, 16> m_array_;
  std::array<std::uint32_t, 16>::iterator m_array_first_;

  static const std::array<std::uint32_t, 64> k_array_;
  static const std::array<std::uint32_t, 64> s_array_;

private:
  static std::uint32_t left_rotate(std::uint32_t x, std::uint32_t c) {
    return (x << c) | (x >> (32 - c));
  }

  template <class OutputIterator>
  static void uint32_to_byte(std::uint32_t n, OutputIterator &first) {

    *first++ = n & 0xff;
    *first++ = (n >> 8) & 0xff;
    *first++ = (n >> 16) & 0xff;
    *first++ = (n >> 24) & 0xff;
  }

  template <class OutputIterator>
  static void uint32_to_hex(std::uint32_t n, OutputIterator &first) {
    const char *hex_chars = "0123456789abcdef";

    std::uint32_t b;

    b = n & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];

    b = (n >> 8) & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];

    b = (n >> 16) & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];

    b = (n >> 24) & 0xff;
    *first++ = hex_chars[b >> 4];
    *first++ = hex_chars[b & 0xf];
  }

private:
  void reset_m_array() { m_array_first_ = m_array_.begin(); }

  template <class InputIterator>
  void bytes_to_m_array(InputIterator &first,
                        std::array<std::uint32_t, 16>::iterator m_array_last) {
    for (; m_array_first_ != m_array_last; ++m_array_first_) {
      *m_array_first_ = *first++;
      *m_array_first_ |= *first++ << 8;
      *m_array_first_ |= *first++ << 16;
      *m_array_first_ |= *first++ << 24;
    }
  }

  template <class InputIterator>
  void true_bit_to_m_array(InputIterator &first, std::ptrdiff_t chunk_length) {
    switch (chunk_length % 4) {
    case 0:
      *m_array_first_++ = 0x00000080;
      break;
    case 1:
      *m_array_first_++ = *first++;
      *m_array_first_ |= 0x00008000;
      break;
    case 2:
      *m_array_first_++ = *first++;
      *m_array_first_ |= *first++ << 8;
      *m_array_first_ |= 0x00800000;
      break;
    case 3:
      *m_array_first_++ = *first++;
      *m_array_first_ |= *first++ << 8;
      *m_array_first_ |= *first++ << 16;
      *m_array_first_ |= 0x80000000;
      break;
    }
  }

  void zeros_to_m_array(std::array<std::uint32_t, 16>::iterator m_array_last) {
    for (; m_array_first_ != m_array_last; ++m_array_first_) {
      *m_array_first_ = 0;
    }
  }

  void original_length_bits_to_m_array(std::uint64_t original_length_bits) {
    original_length_bits &= 0xffffffffffffffff;
    *m_array_first_++ = (original_length_bits) & 0x00000000ffffffff;
    *m_array_first_++ = (original_length_bits & 0xffffffff00000000) >> 32;
  }

  void hash_chunk() {
    std::uint32_t A = a0_;
    std::uint32_t B = b0_;
    std::uint32_t C = c0_;
    std::uint32_t D = d0_;

    std::uint32_t F;
    unsigned int g;

    for (unsigned int i = 0; i < 64; ++i) {
      if (i < 16) {
        F = (B & C) | ((~B) & D);
        g = i;
      } else if (i < 32) {
        F = (D & B) | ((~D) & C);
        g = (5 * i + 1) & 0xf;
      } else if (i < 48) {
        F = B ^ C ^ D;
        g = (3 * i + 5) & 0xf;
      } else {
        F = C ^ (B | (~D));
        g = (7 * i) & 0xf;
      }

      std::uint32_t D_temp = D;
      D = C;
      C = B;
      B += left_rotate(A + F + k_array_[i] + m_array_[g], s_array_[i]);
      A = D_temp;
    }

    a0_ += A;
    b0_ += B;
    c0_ += C;
    d0_ += D;
  }

public:
  template <class InputIterator>
  void update(InputIterator first, InputIterator last) {

    std::uint64_t original_length_bits = std::distance(first, last) * 8;

    std::ptrdiff_t chunk_length;
    while ((chunk_length = std::distance(first, last)) >= 64) {
      reset_m_array();
      bytes_to_m_array(first, m_array_.end());
      hash_chunk();
    }

    reset_m_array();
    bytes_to_m_array(first, m_array_.begin() + chunk_length / 4);
    true_bit_to_m_array(first, chunk_length);

    if (chunk_length >= 56) {
      zeros_to_m_array(m_array_.end());
      hash_chunk();

      reset_m_array();
      zeros_to_m_array(m_array_.end() - 2);
      original_length_bits_to_m_array(original_length_bits);
      hash_chunk();
    } else {
      zeros_to_m_array(m_array_.end() - 2);
      original_length_bits_to_m_array(original_length_bits);
      hash_chunk();
    }
  }

public:
  md5() : a0_(0x67452301), b0_(0xefcdab89), c0_(0x98badcfe), d0_(0x10325476) {}

  template <class Container> void digest(Container &container) {
    container.resize(16);
    auto it = container.begin();

    uint32_to_byte(a0_, it);
    uint32_to_byte(b0_, it);
    uint32_to_byte(c0_, it);
    uint32_to_byte(d0_, it);
  }

  template <class Container> void hex_digest(Container &container) {
    container.resize(32);
    auto it = container.begin();

    uint32_to_hex(a0_, it);
    uint32_to_hex(b0_, it);
    uint32_to_hex(c0_, it);
    uint32_to_hex(d0_, it);
  }
};

const std::array<std::uint32_t, 64> md5::k_array_ = {
    0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a,
    0xa8304613, 0xfd469501, 0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
    0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821, 0xf61e2562, 0xc040b340,
    0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
    0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8,
    0x676f02d9, 0x8d2a4c8a, 0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
    0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70, 0x289b7ec6, 0xeaa127fa,
    0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
    0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92,
    0xffeff47d, 0x85845dd1, 0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
    0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};

const std::array<std::uint32_t, 64> md5::s_array_ = {
    7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
    5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20, 5, 9,  14, 20,
    4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
    6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};

std::string MD5_hash_string(std::string const &data) {
  std::string data_hex_digest;
  md5 hash;
  hash.update(data.begin(), data.end());
  hash.hex_digest(data_hex_digest);
  //
  std::stringstream s;
  s << data_hex_digest;
  std::string converted(s.str());
  return converted;
}

template <typename T> T ConvertHex_to_T(std::string const &data) {
  // Be careful, the T could overflow if T is an "int" though not for
  // a size_t.
  std::vector<std::string> LChar{"0", "1", "2", "3", "4", "5", "6", "7",
                                 "8", "9", "a", "b", "c", "d", "e", "f"};
  auto GetPosition = [&](std::string const &eChar) -> int {
    for (int iPos = 0; iPos < 16; iPos++)
      if (LChar[iPos] == eChar)
        return iPos;
    std::cerr << "Wrong character\n";
    throw TerminalException{1};
  };
  T sum = 0;
  T epow = 1;
  for (size_t u = 0; u < data.size(); u++) {
    int pos = GetPosition(data.substr(u, 1));
    sum += pos * epow;
    epow *= 16;
  }
  return sum;
}

template <typename T> T MD5_hash_T(std::string const &data) {
  std::string data_out = MD5_hash_string(data);
  return ConvertHex_to_T<T>(data_out);
}

// Murmurhash function
// Code from wikipedia

static inline uint32_t murmur_32_scramble(uint32_t k) {
  k *= 0xcc9e2d51;
  k = (k << 15) | (k >> 17);
  k *= 0x1b873593;
  return k;
}

uint32_t murmur3_32(const uint8_t *key, size_t len, uint32_t seed) {
  uint32_t h = seed;
  uint32_t k;
  /* Read in groups of 4. */
  for (size_t i = len >> 2; i; i--) {
    // Here is a source of differing results across endiannesses.
    // A swap here has no effects on hash properties though.
    memcpy(&k, key, sizeof(uint32_t));
    key += sizeof(uint32_t);
    h ^= murmur_32_scramble(k);
    h = (h << 13) | (h >> 19);
    h = h * 5 + 0xe6546b64;
  }
  /* Read the rest. */
  k = 0;
  for (size_t i = len & 3; i; i--) {
    k <<= 8;
    k |= key[i - 1];
  }
  // A swap is *not* necessary here because the preceding loop already
  // places the low bytes in the low places according to whatever endianness
  // we use. Swaps only apply when the memory is copied in a chunk.
  h ^= murmur_32_scramble(k);
  /* Finalize. */
  h ^= uint32_t(len);
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

template <typename T> inline T unaligned_load(void const *ptr) noexcept {
  // using memcpy so we don't get into unaligned load problems.
  // compiler should optimize this very well anyways.
  T t;
  std::memcpy(&t, ptr, sizeof(T));
  return t;
}

// Code adapted from robin_hood hash function
inline size_t robin_hood_hash_bytes(void const *ptr, size_t len,
                                    const uint64_t &seed) noexcept {
  static constexpr uint64_t m = UINT64_C(0xc6a4a7935bd1e995);
  //    static constexpr uint64_t seed = UINT64_C(0xe17a1465);
  static constexpr unsigned int r = 47;

  auto const *const data64 = static_cast<uint64_t const *>(ptr);
  uint64_t h = seed ^ (len * m);

  size_t const n_blocks = len / 8;
  for (size_t i = 0; i < n_blocks; ++i) {
    auto k = unaligned_load<uint64_t>(data64 + i);

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }
  auto const *const data8 =
      reinterpret_cast<uint8_t const *>(data64 + n_blocks);
  switch (len & 7U) {
  case 7:
    h ^= static_cast<uint64_t>(data8[6]) << 48U;
    [[fallthrough]];
  case 6:
    h ^= static_cast<uint64_t>(data8[5]) << 40U;
    [[fallthrough]];
  case 5:
    h ^= static_cast<uint64_t>(data8[4]) << 32U;
    [[fallthrough]];
  case 4:
    h ^= static_cast<uint64_t>(data8[3]) << 24U;
    [[fallthrough]];
  case 3:
    h ^= static_cast<uint64_t>(data8[2]) << 16U;
    [[fallthrough]];
  case 2:
    h ^= static_cast<uint64_t>(data8[1]) << 8U;
    [[fallthrough]];
  case 1:
    h ^= static_cast<uint64_t>(data8[0]);
    h *= m;
    [[fallthrough]];
  default:
    break;
  }

  h ^= h >> r;

  // not doing the final step here, because this will be done by keyToIdx
  // anyways h *= m; h ^= h >> r;
  return static_cast<size_t>(h);
}

template <typename T> std::size_t hash_from_stream(T const &val) {
  std::stringstream s;
  s << val;
  std::string converted(s.str());
  return std::hash<std::string>()(converted);
}

template <typename T>
size_t ComputeHashTwoMap(size_t const &seed,
                         std::map<T, size_t> const &ListDiagNorm,
                         std::map<T, size_t> const &ListOffDiagNorm) {
  auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
    seed ^= new_hash + 0x9e3779b8 + (seed << 6) + (seed >> 2);
  };
  size_t hash = seed;
  auto update_from_map = [&](std::map<T, size_t> const &map) -> void {
    for (auto &kv : map) {
      size_t hash1 = std::hash<T>()(kv.first);
      size_t hash2 = kv.second;
      combine_hash(hash, hash1);
      combine_hash(hash, hash2);
    }
  };
  update_from_map(ListDiagNorm);
  update_from_map(ListOffDiagNorm);
  return hash;
}

namespace std {
template <typename T> struct hash<std::vector<T>> {
  std::size_t operator()(const std::vector<T> &V) const {
    if constexpr (!std::is_arithmetic<T>::value) {
      auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
        seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      size_t seed = 0;
      for (size_t i = 0; i < V.size(); i++) {
        size_t e_hash = std::hash<T>()(V[i]);
        combine_hash(seed, e_hash);
      }
      return seed;
    }
    if constexpr (std::is_arithmetic<T>::value) {
      const T *ptr_T = V.data();
      const uint8_t *ptr_i = reinterpret_cast<const uint8_t *>(ptr_T);
      size_t len = sizeof(T) * V.size();
      uint32_t seed = 0x1b873540;
      return murmur3_32(ptr_i, len, seed);
    }
  }
};
template <typename T1, typename T2> struct hash<std::pair<T1, T2>> {
  std::size_t operator()(const std::pair<T1, T2> &ePair) const {
    auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    size_t seed = std::hash<T1>()(ePair.first);
    size_t e_hash = std::hash<T2>()(ePair.second);
    combine_hash(seed, e_hash);
    return seed;
  }
};
// clang-format off
}  // namespace std
// clang-format on

// clang-format off
#endif  // SRC_BASIC_HASH_FUNCTIONS_H_
// clang-format on
