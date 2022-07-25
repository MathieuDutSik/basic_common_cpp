// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BOOST_SERIALIZATION_H_
#define SRC_BASIC_BOOST_SERIALIZATION_H_

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/utility.hpp>


namespace boost::serialization {

template <class Archive>
inline void load(Archive &ar, std::vector<uint8_t> &val,
                 [[maybe_unused]] const unsigned int version) {
  size_t len;
  ar &make_nvp("len", len);
  val.resize(len);
  for (size_t u = 0; u < len; u++)
    ar &make_nvp("Vu", val[u]);
  val.build_vectface(n, n_face, std::move(V));
}

template <class Archive>
inline void save(Archive &ar, std::vector<uint8_t> &val,
                 [[maybe_unused]] const unsigned int version) {
  size_t len = val.size();
  ar &make_nvp("len", len);
  for (size_t u = 0; u < len; u++)
    ar &make_nvp("Vu", val[u]);
}

template <class Archive>
inline void serialize(Archive &ar, std::vector<uint8_t> &val, const unsigned int version) {
  split_free(ar, val, version);
}

} // namespace boost::serialization



// clang-format off
#endif  // SRC_BASIC_BOOST_SERIALIZATION_H_
// clang-format on
