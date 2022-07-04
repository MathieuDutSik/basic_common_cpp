// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_STRING_H_
#define SRC_BASIC_BASIC_STRING_H_

#include "Temp_common.h"
#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

std::string STRING_GETENV(std::string const &eStr) {
  //  std::string ePrefixAlti;
  //  std::cerr << "STRING_GETENV, step 1\n";
  //  char eStrEnvVar[]="ALTIMETER_DIRECTORY";
  //  std::cerr << "STRING_GETENV, step 2\n";
  char *ePre = std::getenv(eStr.c_str());
  //  std::cerr << "STRING_GETENV, step 3\n";
  if (ePre == NULL) {
    std::cerr << "Error in reading the environment variable : " << eStr << "\n";
    throw TerminalException{1};
  }
  //  std::cerr << "ePre=" << (void*)ePre << "\n";
  std::string eStrRet = ePre;
  //  std::cerr << "STRING_GETENV, step 4\n";
  return eStrRet;
}

bool STRING_IsStringReduceToSpace(std::string const &eStr) {
  size_t len = eStr.length();
  std::string eChar = " ";
  for (size_t i = 0; i < len; i++) {
    std::string eSubChar = eStr.substr(i, 1);
    if (eSubChar != eChar)
      return false;
  }
  return true;
}

int STRING_GetCharPositionInString(std::string const &eStr,
                                   std::string const &eChar) {
  size_t len = eStr.length();
  for (size_t i = 0; i < len; i++) {
    std::string eSubChar = eStr.substr(i, 1);
    if (eSubChar == eChar)
      return static_cast<int>(i);
  }
  return -1;
}

std::string StringSubstitution(
    std::string const &FileIN,
    std::vector<std::pair<std::string, std::string>> const &ListSubst) {
  std::string retStr;
  size_t len = FileIN.size();
  size_t pos = 0;
  std::string CharDollar = "$";
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto GetIPair = [&](size_t const &pos_inp) -> size_t {
    std::vector<size_t> ListMatch;
    for (size_t iPair = 0; iPair < ListSubst.size(); iPair++) {
      std::string eStr = ListSubst[iPair].first;
      size_t lenStr = eStr.size();
      if (pos_inp + lenStr < len) {
        std::string eSubStr = FileIN.substr(pos_inp, lenStr);
        if (eSubStr == eStr)
          ListMatch.push_back(iPair);
      }
    }
    if (ListMatch.size() > 1) {
      std::cerr << "Several strings are matching. Bad situation\n";
      throw TerminalException{1};
    }
    if (ListMatch.size() == 0)
      return miss_val;
    if (ListMatch.size() == 1)
      return ListMatch[0];
    return miss_val;
  };
  while (true) {
    if (pos == len)
      break;
    std::string eChar = FileIN.substr(pos, 1);
    if (eChar == CharDollar) {
      pos++;
      size_t iPair = GetIPair(pos);
      if (iPair == miss_val) {
        retStr += eChar;
      } else {
        retStr += ListSubst[size_t(iPair)].second;
        pos += ListSubst[size_t(iPair)].first.size();
      }
    } else {
      retStr += eChar;
      pos++;
    }
  }
  return retStr;
}

bool IsFullyNumeric(std::string const &eStr) {
  std::string eLS = " 0123456789.";
  size_t nbChar = eStr.size();
  for (size_t uChar = 0; uChar < nbChar; uChar++) {
    std::string eChar = eStr.substr(uChar, 1);
    auto TheFind = [&](std::string const &eCharIn) -> bool {
      size_t nbCase = eLS.size();
      for (size_t iCase = 0; iCase < nbCase; iCase++) {
        std::string eCase = eLS.substr(iCase, 1);
        if (eCase == eCharIn)
          return true;
      }
      return false;
    };
    if (TheFind(eChar) == false)
      return false;
  }
  return true;
}

std::string DoubleTo4dot2f(double const &x) {
  std::stringstream s;
  s << std::fixed;
  s << std::setprecision(2);
  s << x;
  std::string converted(s.str());
  return converted;
}

std::string DoubleTo4dot1f(double const &x) {
  std::stringstream s;
  s << std::fixed;
  s << std::setprecision(1);
  s << x;
  std::string converted(s.str());
  return converted;
}

std::string DoubleToString(double const &x) {
  std::stringstream s;
  s << std::fixed;
  s << std::setprecision(9);
  s << x;
  std::string converted(s.str());
  return converted;
}

std::string IntToString(int const &x) {
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}

int StringToInt(std::string const &str) {
  int val;
  std::istringstream(str) >> val;
  return val;
}

std::string StringVectorStringGAP(std::vector<std::string> const &LStr) {
  std::string ret = "[";
  size_t len = LStr.size();
  for (size_t i = 0; i < len; i++) {
    if (i > 0)
      ret += ",";
    ret += LStr[i];
  }
  ret += "]";
  return ret;
}

std::string LongToString(long const &x) {
  std::stringstream s;
  s << x;
  std::string converted(s.str());
  return converted;
}

int GetNumberDigit(int const &eVal) {
  int nbDigit = 1;
  while (true) {
    if (eVal < pow(10, nbDigit))
      return nbDigit;
    nbDigit++;
  }
}

std::string StringNumber(int const &nb, int const &nbDigit) {
  if (nb < 0) {
    std::cerr << "We have nb < 0 which is an error\n";
    std::cerr << "nb=" << nb << "\n";
    throw TerminalException{1};
  }
  if (nb > pow(10, nbDigit)) {
    std::stringstream s;
    s << "Critical error in StringNumber\n";
    s << "nb=" << nb << "\n";
    s << "nbDigit=" << nbDigit << "\n";
    std::string eStr(s.str());
    throw eStr;
  }
  int idx = 1;
  while (true) {
    if (nb < pow(10, idx)) {
      std::string TheStr;
      for (int i = 0; i < nbDigit - idx; i++)
        TheStr += "0";
      TheStr += std::to_string(nb);
      return TheStr;
    }
    idx++;
  }
}

std::string UpperCaseToLowerCase(std::string const &dataIn) {
  std::string dataRet = dataIn;
  std::transform(dataRet.begin(), dataRet.end(), dataRet.begin(), ::tolower);
  return dataRet;
}

std::string STRING_RemoveSpacesBeginningEnd(std::string const &eStr) {
  size_t len = eStr.size();
  std::vector<int> ListIsSpace(len, 0);
  std::string eSpace = " ";
  for (size_t i = 0; i < len; i++) {
    std::string eChar = eStr.substr(i, 1);
    if (eChar == eSpace)
      ListIsSpace[i] = 1;
  }
  size_t miss_val = std::numeric_limits<size_t>::max();
  size_t PosLow = miss_val;
  for (size_t i = 0; i < len; i++)
    if (PosLow == miss_val)
      if (ListIsSpace[i] == 0)
        PosLow = i;
  size_t PosUpp = miss_val;
  for (size_t i = 0; i < len; i++) {
    size_t j = len - 1 - i;
    if (PosUpp == miss_val)
      if (ListIsSpace[j] == 0)
        PosUpp = j;
  }
  std::string RetStr;
  if (PosLow == miss_val) {
    return RetStr;
  }
  for (size_t iPos = PosLow; iPos <= size_t(PosUpp); iPos++)
    RetStr += eStr.at(iPos);
  return RetStr;
}

template <typename F>
void STRING_Split_f(std::string const &eStrA, std::string const &eStrB, F &f) {
  size_t lenA = eStrA.length();
  size_t lenB = eStrB.length();
  //  std::cerr << "lenA=" << lenA << " lenB=" << lenB << "\n";
  std::vector<int> ListStatus(lenA, 1);
  if (lenA >= lenB) {
    for (size_t iA = 0; iA <= lenA - lenB; iA++)
      if (ListStatus[iA] == 1) {
        auto test = [&]() -> bool {
          for (size_t iB = 0; iB < lenB; iB++)
            if (eStrB[iB] != eStrA[iA + iB])
              return false;
          return true;
        };
        if (test()) {
          for (size_t iB = 0; iB < lenB; iB++) {
            ListStatus[iA + iB] = 0;
          }
        }
      }
  }
  size_t prev_idx = 0;
  size_t iA = 0;
  while (true) {
    if (ListStatus[iA] == 0) {
      if (iA > prev_idx) {
        f(eStrA.substr(prev_idx, iA - prev_idx));
      }
      iA += lenB;
      prev_idx = iA;
    } else {
      iA++;
    }
    if (iA == lenA)
      break;
  }
  if (iA > prev_idx) {
    f(eStrA.substr(prev_idx, iA - prev_idx));
  }
}

// Example of use eStrA="A B C D" and eStrB=" "
std::vector<std::string> STRING_Split(std::string const &eStrA,
                                      std::string const &eStrB) {
  std::vector<std::string> RetList;
  auto f = [&](std::string const &eStr) -> void { RetList.push_back(eStr); };
  STRING_Split_f(eStrA, eStrB, f);
  return RetList;
}

// The String is supposed to be "str0" + hs0 + "str1" + hs1 + "hs2"
// and we return a standard vector of [hs0, hs1]
std::optional<std::vector<std::string>>
STRING_ParseSingleLine_opt(std::string const &strin,
                       std::vector<std::string> const &LStr) {
  size_t len1 = LStr[0].size();
  size_t lentot = strin.size();
  std::string estr = strin.substr(0, len1);
  if (estr != LStr[0]) {
    return {};
  }
  size_t pos = len1;
  size_t nbBlock = LStr.size() - 1;
  std::vector<std::string> LRet(nbBlock);
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto GetInit = [&](std::string const &strSearch,
                     size_t const &posStart) -> size_t {
    size_t lenSearch = strSearch.size();
    for (size_t posi = posStart; posi <= lentot - lenSearch; posi++) {
      std::string strO = strin.substr(posi, lenSearch);
      if (strO == strSearch) {
        return posi;
      }
    }
    return miss_val;
  };
  for (size_t iBlock = 0; iBlock < nbBlock; iBlock++) {
    std::string strSearch = LStr[iBlock + 1];
    size_t lenS = strSearch.size();
    size_t posF = GetInit(strSearch, pos);
    if (posF == miss_val) {
      return {};
    }
    size_t lenO = posF - pos;
    std::string strO = strin.substr(pos, lenO);
    LRet[iBlock] = strO;
    pos += lenO + lenS;
  }
  return LRet;
}

std::vector<int> STRING_Split_Int(std::string const &eStrA,
                                  std::string const &eStrB) {
  std::vector<std::string> LStr = STRING_Split(eStrA, eStrB);
  std::vector<int> LInt;
  for (auto &eStr : LStr) {
    int eVal;
    std::istringstream(eStr) >> eVal;
    LInt.push_back(eVal);
  }
  return LInt;
}

// Example of use eStrA="A B C D" and eStrB=" "
// Difference with the above is that splitting ";;" gives you a list
// of entries as {"", "", ""}, i.e. last and first entry and entry in the middle
std::vector<std::string> STRING_Split_Strict(std::string const &eStrA,
                                             std::string const &eStrB) {
  size_t lenA = eStrA.length();
  size_t lenB = eStrB.length();
  std::vector<size_t> ListStatus(lenA, 0);
  size_t idx = 0;
  for (size_t iA = 0; iA < lenA - lenB; iA++) {
    size_t sumEnt = 0;
    for (size_t iB = 0; iB < lenB; iB++)
      sumEnt += ListStatus[iA + iB];
    if (sumEnt == 0) {
      bool IsMatch = true;
      for (size_t iB = 0; iB < lenB; iB++) {
        std::string eCharA = eStrA.substr(iA + iB, 1);
        std::string eCharB = eStrB.substr(iB, 1);
        if (eCharA != eCharB)
          IsMatch = false;
      }
      if (IsMatch) {
        idx++;
        for (size_t iB = 0; iB < lenB; iB++)
          ListStatus[iA + iB] = idx;
      }
    }
  }
  size_t nbEnt = idx + 1;
  std::vector<std::string> RetList(nbEnt);
  size_t miss_val = std::numeric_limits<size_t>::max();
  for (size_t iEnt = 0; iEnt < nbEnt; iEnt++) {
    size_t posFirst = miss_val, posLast = miss_val;
    if (iEnt == 0) {
      posFirst = 0;
    } else {
      bool WeFound = false;
      for (size_t iA = 0; iA < lenA; iA++) {
        if (!WeFound) {
          if (ListStatus[iA] == iEnt) {
            WeFound = true;
            posFirst = iA + lenB;
          }
        }
      }
    }
    if (iEnt == nbEnt - 1) {
      posLast = lenA - 1;
    } else {
      bool WeFound = false;
      for (size_t iA = 0; iA < lenA; iA++) {
        if (!WeFound) {
          if (ListStatus[iA] == iEnt + 1) {
            WeFound = true;
            posLast = iA - 1;
          }
        }
      }
    }
    if (posFirst == miss_val || posLast == miss_val) {
      std::cerr << "posFirst = " << posFirst << "  posLast = " << posLast
                << "\n";
      std::cerr << "Positions have not been found\n";
      throw TerminalException{1};
    }
    RetList[iEnt] += eStrA.substr(posFirst, 1 + posLast - posFirst);
  }
  return RetList;
}

std::string STRING_Replace(std::string const &eStrA, std::string const &eStrB,
                           std::string const &eStrC) {
  std::vector<std::string> LStr = STRING_Split(eStrA, eStrB);
  std::string str = LStr[0];
  size_t len = LStr.size();
  for (size_t i = 1; i < len; i++)
    str += eStrC + LStr[i];
  return str;
}

std::vector<std::string> STRING_SplitCharNb(std::string const &str) {
  auto IsNumber = [](std::string const &eChar) {
    std::vector<std::string> ListCharNb{"-", "0", "1", "2", "3", "4",
                                        "5", "6", "7", "8", "9"};
    for (auto &fCharNb : ListCharNb)
      if (fCharNb == eChar)
        return true;
    return false;
  };
  size_t len = str.size();
  std::vector<bool> ListStat(len);
  for (size_t i = 0; i < len; i++) {
    std::string eChar = str.substr(i, 1);
    ListStat[i] = IsNumber(eChar);
  }
  std::string TotStr;
  for (size_t i = 0; i < len; i++) {
    TotStr += str.substr(i, 1);
    if (i < len - 1) {
      if (ListStat[i] != ListStat[i + 1]) {
        TotStr += "WRK";
      }
    }
  }
  //  std::cerr << "TotStr=" << TotStr << "\n";
  return STRING_Split(TotStr, "WRK");
}

std::string FILE_GetExtension(std::string const &eFile) {
  std::vector<std::string> LStr = STRING_Split(eFile, "/");
  std::string eFinal = LStr[LStr.size() - 1];
  std::vector<std::string> LBlck = STRING_Split(eFile, ".");
  return LBlck[LBlck.size() - 1];
}

// clang-format off
#endif  // SRC_BASIC_BASIC_STRING_H_
// clang-format on
