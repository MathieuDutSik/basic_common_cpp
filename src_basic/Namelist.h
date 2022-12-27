// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_NAMELIST_H_
#define SRC_BASIC_NAMELIST_H_

#include "Basic_file.h"
#include "Basic_string.h"
#include "Temp_common.h"
#include "hash_functions.h"
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <optional>


struct NamelistException {
  int val;
};

bool NAMELIST_ReadBoolValue(std::string const &eVarValue) {
  if (eVarValue == ".F.")
    return false;
  if (eVarValue == ".T.")
    return true;
  if (eVarValue == "F")
    return false;
  if (eVarValue == "T")
    return true;
  std::cerr << "Boolean value has not been found\n";
  std::cerr << "eVarValue = " << eVarValue << "\n";
  std::cerr << "Allowed: T / F / .T. / .F.\n";
  throw NamelistException{1};
}

std::string NAMELIST_ConvertFortranStringToCppString(std::string const &eStr) {
  int len = eStr.length();
  if (len == 0)
    return "";
  std::string eFirstChar = eStr.substr(0, 1);
  std::string eLastChar = eStr.substr(len - 1, 1);
  int RemovableEnding = 0;
  if (eFirstChar == "'" || eFirstChar == "\"") {
    RemovableEnding = 1;
    if (eFirstChar != eLastChar) {
      std::cerr << "eFirstChar = " << eFirstChar << "\n";
      std::cerr << " eLastChar = " << eLastChar << "\n";
      std::cerr << "The character used for noting beginning and end of string "
                   "should be identical\n";
      throw NamelistException{1};
    }
  }
  if (RemovableEnding == 1) {
    for (int i = 1; i < len - 1; i++) {
      std::string eChar = eStr.substr(i, 1);
      if (eChar == eFirstChar) {
        std::cerr << "It is illegal to have \" or ' character in the middle\n";
        throw TerminalException{1};
      }
    }
    return eStr.substr(1, len - 2);
  }
  return eStr;
}

std::vector<std::string>
NAMELIST_ConvertFortranListStringToCppListString(std::string const &eStr) {
  int len = eStr.length();
  if (len == 0)
    return {};
  std::string eFirstChar = eStr.substr(0, 1);
  std::string eLastChar = eStr.substr(len - 1, 1);
  if (eFirstChar != "'" && eFirstChar != "\"") {
    std::cerr << "eStr=" << eStr << "\n";
    std::cerr
        << "For list of strings, one should use string \"  \"   or '    '   \n";
    throw NamelistException{1};
  }
  if (eLastChar != "'" && eLastChar != "\"") {
    std::cerr << "eStr=" << eStr << "\n";
    std::cerr
        << "For list of strings, one should use string \"  \"   or '    '   \n";
    throw NamelistException{1};
  }
  if (eFirstChar != eLastChar) {
    std::cerr << "eStr=" << eStr << "\n";
    std::cerr << "eFirstChar=" << eFirstChar << "\n";
    std::cerr << "eLastChar=" << eLastChar << "\n";
    std::cerr << "No coherency in endings\n";
    throw NamelistException{1};
  }
  std::string eSepChar = eFirstChar;
  int IsInString = 0;
  std::string eFound = "";
  std::vector<std::string> eListStr;
  for (int i = 0; i < len; i++) {
    std::string eChar = eStr.substr(i, 1);
    if (eChar == eSepChar) {
      eFound += eChar;
      if (IsInString == 1) {
        IsInString = 0;
        std::string eCppStr = NAMELIST_ConvertFortranStringToCppString(eFound);
        eListStr.push_back(eCppStr);
        eFound = "";
      } else {
        IsInString = 1;
      }
    } else {
      if (IsInString == 1)
        eFound += eChar;
    }
  }
  return eListStr;
}

std::vector<double> NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(
    std::string const &eVarValue) {
  std::string eSepChar = ",";
  std::vector<std::string> ListStr = STRING_Split(eVarValue, eSepChar);
  std::vector<double> eListRetDouble;
  int siz = ListStr.size();
  for (int i = 0; i < siz; i++) {
    std::string eStr1 = ListStr[i];
    std::string eStr2 = STRING_RemoveSpacesBeginningEnd(eStr1);
    double eVal = ParseScalar<double>(eStr2);
    eListRetDouble.push_back(eVal);
  }
  return eListRetDouble;
}

std::vector<int> NAMELIST_ConvertFortranStringListIntToCppVectorInt(
    std::string const &eVarValue) {
  std::string eSepChar = ",";
  std::vector<std::string> ListStr = STRING_Split(eVarValue, eSepChar);
  std::vector<int> eListRetInt;
  int siz = ListStr.size();
  for (int i = 0; i < siz; i++) {
    std::string eStr1 = ListStr[i];
    std::string eStr2 = STRING_RemoveSpacesBeginningEnd(eStr1);
    int eVal = ParseScalar<int>(eStr2);
    eListRetInt.push_back(eVal);
  }
  return eListRetInt;
}

std::optional<std::string> get_default(std::string const& strin) {
  std::string prefix = "Default: ";
  size_t prefix_s = prefix.size();
  if (strin.size() < prefix_s) {
    return {};
  }
  if (strin.substr(0, prefix_s) != prefix) {
    return {};
  }
  size_t n_char = strin.size();
  for (size_t i_char=0; i_char<n_char; i_char++) {
    std::string e_char = strin.substr(i_char,1);
    if (e_char != "\n") {
      return strin.substr(prefix_s, i_char - prefix_s);
    }
  }
  return strin.substr(prefix_s, n_char - prefix_s);
}

struct SingleBlock {
public:
  std::map<std::string, int> ListIntValues;
  std::map<std::string, bool> ListBoolValues;
  std::map<std::string, double> ListDoubleValues;
  std::map<std::string, std::vector<double>> ListListDoubleValues;
  std::map<std::string, std::vector<int>> ListListIntValues;
  std::map<std::string, std::string> ListStringValues;
  std::map<std::string, std::vector<std::string>> ListListStringValues;
  std::map<std::string, std::string> ListIntValues_doc;
  std::map<std::string, std::string> ListBoolValues_doc;
  std::map<std::string, std::string> ListDoubleValues_doc;
  std::map<std::string, std::string> ListListDoubleValues_doc;
  std::map<std::string, std::string> ListListIntValues_doc;
  std::map<std::string, std::string> ListStringValues_doc;
  std::map<std::string, std::string> ListListStringValues_doc;
  std::vector<std::string> ListNoDefault;
  void setListIntValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListIntValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        ListIntValues[kv.first] = ParseScalar<int>(*opt);
      } else {
        ListIntValues[kv.first] = 0;
        ListNoDefault.push_back(kv.first);
      }
    }
  }
  void setListBoolValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListBoolValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        try {
          ListBoolValues[kv.first] = NAMELIST_ReadBoolValue(*opt);
        }
        catch (NamelistException &e) {
          std::cerr << "Error parsing the boolean kv.second=" << kv.second << " *opt=" << *opt << " |*opt|=" << opt->size() << "\n";
          throw TerminalException{1};
        }
      } else {
        ListBoolValues[kv.first] = false;
        ListNoDefault.push_back(kv.first);
      }
    }
  }
  void setListDoubleValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListDoubleValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        ListDoubleValues[kv.first] = ParseScalar<double>(*opt);
      } else {
        ListDoubleValues[kv.first] = 0;
        ListNoDefault.push_back(kv.first);
      }
    }
  }
  void setListStringValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListStringValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        try {
          ListStringValues[kv.first] = NAMELIST_ConvertFortranStringToCppString(*opt);
        }
        catch (NamelistException &e) {
          std::cerr << "Error parsing the string kv.second=" << kv.second << "\n";
          throw TerminalException{1};
        }
      } else {
        ListStringValues[kv.first] = "";
        ListNoDefault.push_back(kv.first);
      }
    }
  }
  void setListListDoubleValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListListDoubleValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        try {
          ListListDoubleValues[kv.first] = NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(*opt);
        }
        catch (NamelistException &e) {
          std::cerr << "Error parsing the string kv.second=" << kv.second << "\n";
          throw TerminalException{1};
        }
      } else {
        ListListDoubleValues[kv.first] = {};
        ListNoDefault.push_back(kv.first);
      }
    }
  }
  void setListListIntValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListListIntValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        try {
          ListListIntValues[kv.first] = NAMELIST_ConvertFortranStringListIntToCppVectorInt(*opt);
        }
        catch (NamelistException &e) {
          std::cerr << "Error parsing the string kv.second=" << kv.second << "\n";
          throw TerminalException{1};
        }
      } else {
        ListListIntValues[kv.first] = {};
        ListNoDefault.push_back(kv.first);
      }
    }
  }
  void setListListStringValues(std::map<std::string, std::string> const& m) {
    for (auto & kv : m) {
      ListListStringValues_doc[kv.first] = kv.second;
      std::optional<std::string> opt = get_default(kv.second);
      if (opt) {
        try {
          ListListStringValues[kv.first] = NAMELIST_ConvertFortranListStringToCppListString(*opt);
        }
        catch (NamelistException &e) {
          std::cerr << "Error parsing the string kv.second=" << kv.second << "\n";
          throw TerminalException{1};
        }
      } else {
        ListListStringValues[kv.first] = {};
        ListNoDefault.push_back(kv.first);
      }
    }
  }

};


struct FullNamelist {
  std::map<std::string, SingleBlock> ListBlock;
  std::string FileName;
};

std::vector<std::string> ExtractMatchingBool(SingleBlock const &eBlock) {
  std::vector<std::string> ListMatch;
  for (auto & kv : eBlock.ListBoolValues)
    if (kv.second)
      ListMatch.push_back(kv.first);
  return ListMatch;
}

std::string
NAMELIST_FindPositionVariableInBlock(std::string const &FullVarName,
                                     SingleBlock const &eSingleBlock) {
  std::vector<std::string> LStr = STRING_Split(FullVarName, ":");
  std::string eVarName = LStr[0];
  if (eSingleBlock.ListIntValues.count(eVarName) > 0)
    return "int";
  if (eSingleBlock.ListBoolValues.count(eVarName) > 0)
    return "bool";
  if (eSingleBlock.ListDoubleValues.count(eVarName) > 0)
    return "double";
  if (eSingleBlock.ListListDoubleValues.count(eVarName) > 0)
    return "listdouble";
  if (eSingleBlock.ListListIntValues.count(eVarName) > 0)
    return "listint";
  if (eSingleBlock.ListStringValues.count(eVarName) > 0)
    return "string";
  if (eSingleBlock.ListListStringValues.count(eVarName) > 0)
    return "liststring";
  return "not found";
}

std::string NAMELIST_RemoveAfterCommentChar(std::string const &eStr,
                                            std::string const &eChar) {
  bool WeFound = false;
  std::string RetStr;
  int len = eStr.size();
  for (int i = 0; i < len; i++) {
    std::string fChar = eStr.substr(i, 1);
    if (fChar == eChar)
      WeFound = true;
    if (!WeFound)
      RetStr += eStr.at(i);
  }
  return RetStr;
}

std::string NAMELIST_RemoveAfterLastChar(std::string const &eStr,
                                         std::string const &eLastChar) {
  int iPos = -1;
  int len = eStr.size();
  for (int i = 0; i < len; i++) {
    int j = len - 1 - i;
    if (iPos == -1) {
      std::string eChar = eStr.substr(j, 1);
      if (eChar == eLastChar)
        iPos = j;
    }
  }
  if (iPos == -1)
    return eStr;
  return eStr.substr(0, iPos);
}

std::string NAMELIST_ClearEndOfLine(std::string const &eStr) {
  std::string eCharCommentB = "!";
  std::string eStr3 = NAMELIST_RemoveAfterLastChar(eStr, eCharCommentB);
  //
  int iPos = -1;
  int len = eStr3.size();
  std::string eLastChar = ",";
  for (int i = 0; i < len; i++) {
    int j = len - 1 - i;
    if (iPos == -1) {
      std::string eChar = eStr3.substr(j, 1);
      if (eChar == eLastChar)
        iPos = j;
    }
  }
  if (iPos == -1)
    return eStr3;
  std::string eStrPrior = eStr3.substr(0, iPos);
  std::string eStrPosterior = eStr3.substr(iPos + 1, len - iPos - 1);
  bool test = STRING_IsStringReduceToSpace(eStrPosterior);
  if (test)
    return eStrPrior;
  return eStr3;
}

void NAMELIST_WriteBlock_Kernel(std::ostream &os, std::string const &eBlockName,
                                SingleBlock const &eBlock, bool const& WithDoc) {
  os << "&" << eBlockName << "\n";
  //
  // Integer values
  //
  for (auto & kv : eBlock.ListIntValues) {
    auto iter = eBlock.ListIntValues_doc.find(kv.first);
    if (iter == eBlock.ListIntValues_doc.end() || !WithDoc) {
      os << "  " << kv.first << " = " << kv.second << "\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  //
  // Bool values
  //
  for (auto & kv : eBlock.ListBoolValues) {
    auto iter = eBlock.ListBoolValues_doc.find(kv.first);
    if (iter == eBlock.ListBoolValues_doc.end() || !WithDoc) {
      bool eVal = kv.second;
      std::string eValStr;
      if (!eVal)
        eValStr = "F";
      else
        eValStr = "T";
      os << "  " << kv.first << " = " << eValStr << "\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  //
  // Double values
  //
  for (auto & kv : eBlock.ListDoubleValues) {
    auto iter = eBlock.ListDoubleValues_doc.find(kv.first);
    if (iter == eBlock.ListDoubleValues_doc.end() || !WithDoc) {
      os << "  " << kv.first << " = " << kv.second << "\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  //
  // ListDouble values
  //
  for (auto & kv : eBlock.ListListDoubleValues) {
    auto iter = eBlock.ListListDoubleValues_doc.find(kv.first);
    if (iter == eBlock.ListListDoubleValues_doc.end() || !WithDoc) {
      os << "  " << kv.first << " = ";
      std::vector<double> const& eListDoubl = kv.second;
      int nbDoubl = eListDoubl.size();
      for (int iDoubl = 0; iDoubl < nbDoubl; iDoubl++) {
        if (iDoubl > 0)
          os << ", ";
        os << eListDoubl[iDoubl];
      }
      os << "\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  //
  // ListInt values
  //
  for (auto & kv : eBlock.ListListIntValues) {
    auto iter = eBlock.ListListIntValues_doc.find(kv.first);
    if (iter == eBlock.ListListIntValues_doc.end() || !WithDoc) {
      os << "  " << kv.first << " = ";
      std::vector<int> const& eListInt = kv.second;
      int nbInt = eListInt.size();
      for (int iInt = 0; iInt < nbInt; iInt++) {
        if (iInt > 0)
          os << ", ";
        os << eListInt[iInt];
      }
      os << "\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  //
  // String values
  //
  for (auto & kv : eBlock.ListStringValues) {
    auto iter = eBlock.ListStringValues_doc.find(kv.first);
    if (iter == eBlock.ListStringValues_doc.end() || !WithDoc) {
      os << "  " << kv.first << " = \"" << kv.second << "\"\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  //
  // ListString values
  //
  for (auto & kv : eBlock.ListListStringValues) {
    auto iter = eBlock.ListListStringValues_doc.find(kv.first);
    if (iter == eBlock.ListListStringValues_doc.end() || !WithDoc) {
      os << "  " << kv.first << " = ";
      std::vector<std::string> const& eListStr = kv.second;
      int nbString = eListStr.size();
      for (int iString = 0; iString < nbString; iString++) {
        if (iString > 0)
          os << ", ";
        os << "\"" << eListStr[iString] << "\"";
      }
      os << "\n";
    } else {
      os << "  " << kv.first << " : " << iter->second << "\n";
    }
  }
  os << "/\n";
}

void NAMELIST_WriteBlock(std::ostream &os, std::string const &eBlockName,
                         SingleBlock const &eBlock) {
  NAMELIST_WriteBlock_Kernel(os, eBlockName, eBlock, false);
}


void NAMELIST_WriteNamelistFile(std::ostream &os,
                                FullNamelist const &eFull) {
  int iBlock = 0;
  for (auto & kv : eFull.ListBlock) {
    std::string const& eBlockName = kv.first;
    SingleBlock const& eBlock = kv.second;
    if (iBlock > 0)
      os << "\n\n";
    NAMELIST_WriteBlock(os, eBlockName, eBlock);
    iBlock++;
  }
}

std::vector<std::string>
NAMELIST_ListTrueEntryBool(FullNamelist const &eFull,
                           std::string const &eBlockName) {
  std::vector<std::string> ListString;
  for (auto &kv : eFull.ListBlock.at(eBlockName).ListBoolValues)
    if (kv.second)
      ListString.push_back(kv.first);
  return ListString;
}

void NAMELIST_ReadNamelistStream(std::istream & is,
                                 FullNamelist &eFull) {
  std::unordered_set<std::pair<std::string, std::string>> ListInsertValues;
  auto parsing_error_end = [&](std::string const &eBlockName,
                               std::string const &eVarName,
                               std::string const &TypeVar) -> void {
    std::cerr << "Error reading in the block " << eBlockName << "\n";
    std::cerr << "Variable eVarName=" << eVarName << " should be a " << TypeVar
              << "\n";
    std::cerr << "Please correct your input file\n";
    throw TerminalException{1};
  };
  bool InBlock = false;
  std::string eBlockName;
  std::map<std::string, std::set<std::string>> ls_string;
  while (!is.eof()) {
    std::string Ampersand = "&";
    std::string strTab = "\t";
    std::string PreStr;
    std::getline(is, PreStr);
    std::string eCharComment = "!";
    std::string PreStrB = NAMELIST_RemoveAfterCommentChar(PreStr, eCharComment);
    std::string eStr = STRING_RemoveSpacesBeginningEnd(PreStrB);
    int len = eStr.length();
    if (eStr.find(strTab) != std::string::npos) {
      std::cerr << "Tabs are not allowed\n";
      std::cerr << "LINE=" << eStr << "\n";
      throw TerminalException{1};
    }
    if (len > 0) {
      if (eStr.find(Ampersand) != std::string::npos) {
        std::string eFirstChar = eStr.substr(0, 1);
        if (eFirstChar != "&") {
          std::cerr << "Error while processing stream\n";
          std::cerr
              << "Error, Ampersand (&) should be only in the first character\n";
          std::cerr << "LINE=" << eStr << "\n";
          throw TerminalException{1};
        }
        std::string strRed = eStr.substr(1, len - 1);
        if (!InBlock) {
          eBlockName = strRed;
          if (eFull.ListBlock.count(eBlockName) == 0) {
            std::cerr << "Find BlockName = " << eBlockName << "\n";
            std::cerr << "which is not in the authorized list\n";
            std::cerr << "LINE=" << eStr << "\n";
            std::cerr << "List of authorized block names:\n";
            for (auto &eBlock : eFull.ListBlock)
              std::cerr << "Block name=" << eBlock.first << "\n";
            throw TerminalException{1};
          }
          InBlock = true;
        } else {
          if (strRed != "END") {
            std::cerr << "Ampersand detected. We should leave with a END\n";
            std::cerr << "LINE=" << eStr << "\n";
            throw TerminalException{1};
          }
          InBlock = false;
        }
      } else {
        if (eStr != "/") {
          std::string eStr3 = NAMELIST_ClearEndOfLine(eStr);
          std::string strEqual = "=";
          int posEqual = STRING_GetCharPositionInString(eStr3, strEqual);
          if (posEqual != -1) {
            int len3 = eStr3.length();
            std::string eStrPrior = eStr3.substr(0, posEqual);
            std::string eStrPosterior =
                eStr3.substr(posEqual + 1, len3 - posEqual - 1);
            std::string eVarName = STRING_RemoveSpacesBeginningEnd(eStrPrior);
            std::pair<std::string, std::string> ePair{eBlockName, eVarName};
            if (ListInsertValues.count(ePair) > 0) {
              std::cerr << "In the block " << eBlockName << "\n";
              std::cerr << "the entry " << eVarName << "\n";
              std::cerr << "is defined two times\n";
              throw TerminalException{1};
            }
            ListInsertValues.insert(ePair);
            //
            std::string eVarValue =
                STRING_RemoveSpacesBeginningEnd(eStrPosterior);
            std::string eVarNature = NAMELIST_FindPositionVariableInBlock(
                eVarName, eFull.ListBlock[eBlockName]);
            if (eVarNature == "not found") {
              NAMELIST_WriteBlock(std::cerr, eBlockName,
                                  eFull.ListBlock[eBlockName]);
              std::cerr << "Error in reading the NAMELIST file. See above "
                           "allowed entries\n";
              std::cerr << "The variable " << eVarName << "\n";
              std::cerr << "is in block " << eBlockName << "\n";
              std::cerr << "but it is not allowed for the chosen application\n";
              throw TerminalException{1};
            }
            if (eVarNature == "int") {
              int eVal = ParseScalar<int>(eVarValue);
              eFull.ListBlock[eBlockName].ListIntValues[eVarName] = eVal;
              ls_string[eBlockName].insert(eVarName);
            }
            if (eVarNature == "bool") {
              try {
                bool eVal = NAMELIST_ReadBoolValue(eVarValue);
                eFull.ListBlock[eBlockName].ListBoolValues[eVarName] =
                    eVal;
                ls_string[eBlockName].insert(eVarName);
              } catch (NamelistException &e) {
                parsing_error_end(eBlockName, eVarName, "bool");
              }
            }
            if (eVarNature == "double") {
              double eVal = ParseScalar<double>(eVarValue);
              eFull.ListBlock[eBlockName].ListDoubleValues[eVarName] = eVal;
              ls_string[eBlockName].insert(eVarName);
            }
            if (eVarNature == "string") {
              try {
                std::string eVal =
                    NAMELIST_ConvertFortranStringToCppString(eVarValue);
                eFull.ListBlock[eBlockName].ListStringValues[eVarName] = eVal;
                ls_string[eBlockName].insert(eVarName);
              } catch (NamelistException &e) {
                parsing_error_end(eBlockName, eVarName, "string");
              }
            }
            if (eVarNature == "listdouble") {
              std::vector<double> eVal =
                  NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(
                      eVarValue);
              eFull.ListBlock[eBlockName]
                  .ListListDoubleValues[eVarName] = eVal;
              ls_string[eBlockName].insert(eVarName);
            }
            if (eVarNature == "listint") {
              std::vector<int> eVal =
                  NAMELIST_ConvertFortranStringListIntToCppVectorInt(eVarValue);
              eFull.ListBlock[eBlockName].ListListIntValues[eVarName] =
                  eVal;
              ls_string[eBlockName].insert(eVarName);
            }
            if (eVarNature == "liststring") {
              try {
                std::vector<std::string> eVal =
                    NAMELIST_ConvertFortranListStringToCppListString(eVarValue);
                eFull.ListBlock[eBlockName]
                    .ListListStringValues[eVarName] = eVal;
                ls_string[eBlockName].insert(eVarName);
              } catch (NamelistException &e) {
                parsing_error_end(eBlockName, eVarName, "liststring");
              }
            }
          } else {
            if (eStr3.size() != 0) {
              std::cerr << "If lines has no = sign then it should be empty\n";
              std::cerr << "str=" << PreStr << "\n";
              throw TerminalException{1};
            }
          }
        } else {
          InBlock = false;
        }
      }
    }
  }
  if (InBlock) {
    std::cerr
        << "Error. When leaving namelist reading, we should be out of block\n";
    throw TerminalException{1};
  }
  for (auto & kv : eFull.ListBlock) {
    auto the_set = ls_string[kv.first];
    for (auto & KeyNotDefault : kv.second.ListNoDefault) {
      auto iter = the_set.find(KeyNotDefault);
      if (iter == the_set.end()) {
        std::cerr << "The key " << KeyNotDefault << " has not been assigned\n";
        std::cerr << "This is needed because it is not a default key\n";
        throw TerminalException{1};
      }
    }
  }
}

void NAMELIST_ReadNamelistFile(std::string const &eFileName,
                               FullNamelist &eFull) {
  if (!IsExistingFile(eFileName)) {
    std::cerr << "The following namelist file is missing\n";
    std::cerr << "eFileName = " << eFileName << "\n";
    throw TerminalException{1};
  }
  std::ifstream INfs(eFileName);
  NAMELIST_ReadNamelistStream(INfs, eFull);
}



void NAMELIST_ReadListString(FullNamelist &eFull, std::vector<std::string> const &ListString) {
  std::string str_tot;
  for (auto const &eStr : ListString) {
    str_tot += eStr;
    str_tot += "\n";
  }
  std::istringstream is(str_tot);
  NAMELIST_ReadNamelistStream(is, eFull);
}



// clang-format off
#endif  // SRC_BASIC_NAMELIST_H_
// clang-format on
