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
#include <vector>

namespace boost::serialization {

  // std::vector

  template <class Archive, typename T>
  inline void serialize(Archive &ar, std::vector<T> &val,
                        [[maybe_unused]] const unsigned int version) {
    size_t len = val.size();
    ar &make_nvp("len", len);
    val.resize(len);
    // always save/load row-major
    for (size_t u = 0; u < len; u++)
      ar &make_nvp("Vu", val[u]);
  }

  // std::unordered_set

  template <class Archive, typename T>
  inline void load(Archive &ar, std::unordered_set<T> &val,
                   [[maybe_unused]] const unsigned int version) {
    size_t n;
    ar &make_nvp("n", n);
    val.clear();
    for (size_t u = 0; u < n; u++) {
      T key;
      ar &make_nvp("key", key);
      val.insert(key);
    }
  }

  template <class Archive, typename T>
  inline void save(Archive &ar, std::unordered_set<T> const &val,
                   [[maybe_unused]] const unsigned int version) {
    size_t n = val.size();
    ar &make_nvp("n", n);
    for (auto& key : val) {
      ar &make_nvp("key", key);
    }
  }

  template <class Archive, typename T>
  inline void serialize(Archive &ar, std::unordered_set<T> &val, const unsigned int version) {
    split_free(ar, val, version);
  }

  // std::unordered_map

  template <class Archive, typename K, typename V>
  inline void load(Archive &ar, std::unordered_map<K, V> &val,
                   [[maybe_unused]] const unsigned int version) {
    size_t n;
    ar &make_nvp("n", n);
    val.clear();
    for (size_t u = 0; u < n; u++) {
      K key;
      V value;
      ar &make_nvp("key", key);
      ar &make_nvp("value", value);
      val.insert(std::make_pair(key, value));
    }
  }

  template <class Archive, typename K, typename V>
  inline void save(Archive &ar, std::unordered_map<K, V> const &val,
                   [[maybe_unused]] const unsigned int version) {
    size_t n = val.size();
    ar &make_nvp("n", n);
    for (auto& kv : val) {
      ar &make_nvp("key", kv.first);
      ar &make_nvp("value", kv.second);
    }
  }

  template <class Archive, typename K, typename V>
  inline void serialize(Archive &ar, std::unordered_map<K, V> &val, const unsigned int version) {
    split_free(ar, val, version);
  }

  // std::optional

  template <class Archive, typename T>
  inline void load(Archive &ar, std::optional<T> &val,
                   [[maybe_unused]] const unsigned int version) {
    bool test;
    ar &make_nvp("has_value", test);
    if (test) {
      T value;
      ar &make_nvp("value", value);
      val = value;
    }
  }

  template <class Archive, typename T>
  inline void save(Archive &ar, std::optional<T> const &val,
                   [[maybe_unused]] const unsigned int version) {
    bool test = val.has_value();
    ar &make_nvp("has_value", test);
    if (test) {
      ar &make_nvp("value", *val);
    }
  }

  template <class Archive, typename T>
  inline void serialize(Archive &ar, std::optional<T> &val, const unsigned int version) {
    split_free(ar, val, version);
  }

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_BASIC_BOOST_SERIALIZATION_H_
// clang-format on
