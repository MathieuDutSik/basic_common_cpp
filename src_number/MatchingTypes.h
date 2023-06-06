// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_MATCHINGTYPES_H_
#define SRC_NUMBER_MATCHINGTYPES_H_

template <typename T, typename Tinf> constexpr T get_signed_max_range() {
  T val_max = std::numeric_limits<T>::max();
  T val_min = std::numeric_limits<T>::min();
  Tinf val_max_tinf = UniversalScalarConversion<Tinf, T>(val_max);
  Tinf val_min_tinf = UniversalScalarConversion<Tinf, T>(val_min);
  if (val_max_tinf < -val_min_tinf) {
    return val_max;
  } else {
    return -val_min;
  }
}

// Find the maximal possible rang to run the program.
// Tinf should be an infinite integer type like mpz_class
template <typename Tinf> std::string get_matching_types(Tinf const &max_val) {
  int16_t val16 = get_signed_max_range<int16_t, Tinf>();
  int32_t val32 = get_signed_max_range<int32_t, Tinf>();
  int64_t val64 = get_signed_max_range<int64_t, Tinf>();
  Tinf val16_T = UniversalScalarConversion<Tinf, int16_t>(val16);
  Tinf val32_T = UniversalScalarConversion<Tinf, int16_t>(val32);
  Tinf val64_T = UniversalScalarConversion<Tinf, int16_t>(val64);
  if (max_val < val16_T)
    return "int16_t";
  if (max_val < val32_T)
    return "int32_t";
  if (max_val < val64_T)
    return "int64_t";
  return "infinite";
}

// clang-format off
#endif  // SRC_NUMBER_MATCHINGTYPES_H_
// clang-format on
