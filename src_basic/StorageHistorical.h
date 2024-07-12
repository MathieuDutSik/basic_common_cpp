// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_STORAGEHISTORICAL_H_
#define SRC_BASIC_STORAGEHISTORICAL_H_

#include <utility>
#include <vector>

template <typename T> struct StorageHistorical {
  struct FullRelInfo {
    bool status;
    uint64_t eTime;
    T eVal;
  };

public:
  // no copy
  StorageHistorical(const StorageSpaceLastN<T> &) = delete;

  // no assign
  StorageHistorical &operator=(const StorageHistorical<T> &) = delete;

  // no move
  StorageHistorical(StorageHistorical<T> &&) = delete;

  StorageHistorical() {
    ThePeriod = -1;
    siz = 0;
  }

  StorageHistorical(int const &_ThePeriod) { ThePeriod = _ThePeriod; }

  void SetPeriod(int const &_ThePeriod) { ThePeriod = _ThePeriod; }

  std::vector<T> RetrieveLastRelevantValues() const {
    std::vector<T> ListVal;
    for (auto &eFull : ListFull)
      if (eFull.status)
        ListVal.push_back(eFull.eVal);
    return ListVal;
  }

  void InsertOneValue(uint64_t const &eTime, T const &eVal) {
    if (ThePeriod < 0) {
      return;
    }
    bool HasFoundPlace = false;
    for (int i = 0; i < siz; i++) {
      int deltaTime = static_cast<int>(eTime - ListFull[i].eTime);
      if (deltaTime > ThePeriod)
        ListFull[i].status = false;
      if (!HasFoundPlace && !ListFull[i].status) {
        ListFull[i] = {true, eTime, eVal};
        HasFoundPlace = true;
      }
    }
    if (!HasFoundPlace) {
      ListFull.push_back({true, eTime, eVal});
      siz++;
    }
  }

private:
  std::vector<FullRelInfo> ListFull;
  int siz;
  int ThePeriod;
};

// clang-format off
#endif  // SRC_BASIC_STORAGEHISTORICAL_H_
// clang-format on
