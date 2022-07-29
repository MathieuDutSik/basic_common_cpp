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

struct SingleBlock {
  std::map<std::string, int> ListIntValues;
  std::map<std::string, bool> ListBoolValues;
  std::map<std::string, double> ListDoubleValues;
  std::map<std::string, std::vector<double>> ListListDoubleValues;
  std::map<std::string, std::vector<int>> ListListIntValues;
  std::map<std::string, std::string> ListStringValues;
  std::map<std::string, std::vector<std::string>> ListListStringValues;
};

struct FullNamelist {
  std::map<std::string, SingleBlock> ListBlock;
  std::string FileName;
};

std::vector<std::string> ExtractMatchingBool(SingleBlock const &eBlock) {
  std::map<std::string, bool>::const_iterator iter =
      eBlock.ListBoolValues.begin();
  std::vector<std::string> ListMatch;
  while (iter != eBlock.ListBoolValues.end()) {
    if (iter->second)
      ListMatch.push_back(iter->first);
    iter++;
  }
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
  std::string RetStr;
  for (int i = 0; i < iPos; i++)
    RetStr += eStr.at(i);
  return RetStr;
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
    double eVal;
    std::istringstream(eStr2) >> eVal;
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
    int eVal;
    std::istringstream(eStr2) >> eVal;
    eListRetInt.push_back(eVal);
  }
  return eListRetInt;
}

void NAMELIST_WriteBlock(std::ostream &os, std::string const &eBlockName,
                         SingleBlock const &eBlock) {
  os << "&" << eBlockName << "\n";
  for (auto & kv : eBlock.ListIntValues)
    os << "  " << kv.first << " = " << kv.second << "\n";
  for (auto & kv : eBlock.ListBoolValues) {
    bool eVal = kv.second;
    std::string eValStr;
    if (!eVal)
      eValStr = "F";
    else
      eValStr = "T";
    os << "  " << kv.first << " = " << eValStr << "\n";
  }
  for (auto & kv : eBlock.ListDoubleValues)
    os << "  " << kv.first << " = " << kv.second << "\n";
  for (auto & kv : eBlock.ListListDoubleValues) {
    os << "  " << kv.first << " = ";
    std::vector<double> const& eListDoubl = kv.second;
    int nbDoubl = eListDoubl.size();
    for (int iDoubl = 0; iDoubl < nbDoubl; iDoubl++) {
      if (iDoubl > 0)
        os << ", ";
      os << eListDoubl[iDoubl];
    }
    os << "\n";
  }
  for (auto & kv : eBlock.ListListIntValues) {
    os << "  " << kv.first << " = ";
    std::vector<int> const& eListInt = kv.second;
    int nbInt = eListInt.size();
    for (int iInt = 0; iInt < nbInt; iInt++) {
      if (iInt > 0)
        os << ", ";
      os << eListInt[iInt];
    }
    os << "\n";
  }
  for (auto & kv : eBlock.ListStringValues)
    os << "  " << kv.first << " = \"" << kv.second << "\"\n";
  for (auto & kv : eBlock.ListListStringValues) {
    os << "  " << kv.first << " = ";
    std::vector<std::string> const& eListStr = kv.second;
    int nbString = eListStr.size();
    for (int iString = 0; iString < nbString; iString++) {
      if (iString > 0)
        os << ", ";
      os << "\"" << eListStr[iString] << "\"";
    }
    os << "\n";
  }
  os << "/\n";
}

void NAMELIST_WriteNamelistFile(std::ostream &os,
                                FullNamelist const &eFullNamelist) {
  int iBlock = 0;
  for (std::map<std::string, SingleBlock>::const_iterator itB =
           eFullNamelist.ListBlock.begin();
       itB != eFullNamelist.ListBlock.end(); ++itB) {
    std::string eBlockName = itB->first;
    SingleBlock eBlock = itB->second;
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

void NAMELIST_ReadNamelistFile(std::string const &eFileName,
                               FullNamelist &eFullNamelist) {
  std::unordered_set<std::pair<std::string, std::string>> ListInsertValues;
  if (!IsExistingFile(eFileName)) {
    std::cerr << "The following namelist file is missing\n";
    std::cerr << "eFileName = " << eFileName << "\n";
    throw TerminalException{1};
  }
  auto parsing_error_end = [&](std::string const &eBlockName,
                               std::string const &eVarName,
                               std::string const &TypeVar) -> void {
    std::cerr << "Error reading in the block " << eBlockName << "\n";
    std::cerr << "Variable eVarName=" << eVarName << " should be a " << TypeVar
              << "\n";
    std::cerr << "Please correct your input file\n";
    throw TerminalException{1};
  };
  std::ifstream INfs(eFileName);
  bool InBlock = false;
  std::string eBlockName;
  while (!INfs.eof()) {
    std::string Ampersand = "&";
    std::string strTab = "\t";
    std::string PreStr;
    std::getline(INfs, PreStr);
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
          std::cerr << "Error while reading namelist file = " << eFileName
                    << "\n";
          std::cerr
              << "Error, Ampersand (&) should be only in the first character\n";
          std::cerr << "LINE=" << eStr << "\n";
          throw TerminalException{1};
        }
        std::string strRed = eStr.substr(1, len - 1);
        if (!InBlock) {
          eBlockName = strRed;
          if (eFullNamelist.ListBlock.count(eBlockName) == 0) {
            std::cerr << "Find BlockName = " << eBlockName << "\n";
            std::cerr << "which is not in the authorized list\n";
            std::cerr << "LINE=" << eStr << "\n";
            std::cerr << "List of authorized block names:\n";
            for (auto &eBlock : eFullNamelist.ListBlock)
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
              std::cerr << "eFileName = " << eFileName << "\n";
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
                eVarName, eFullNamelist.ListBlock[eBlockName]);
            if (eVarNature == "not found") {
              NAMELIST_WriteBlock(std::cerr, eBlockName,
                                  eFullNamelist.ListBlock[eBlockName]);
              std::cerr << "Error in reading the NAMELIST file. See above "
                           "allowed entries\n";
              std::cerr << "The variable " << eVarName << "\n";
              std::cerr << "is in block " << eBlockName << "\n";
              std::cerr << "of the file " << eFileName << "\n";
              std::cerr << "but it is not allowed for the chosen application\n";
              throw TerminalException{1};
            }
            if (eVarNature == "int") {
              int eVal;
              std::istringstream(eVarValue) >> eVal;
              eFullNamelist.ListBlock[eBlockName].ListIntValues[eVarName] =
                  eVal;
            }
            if (eVarNature == "bool") {
              try {
                bool eVal = NAMELIST_ReadBoolValue(eVarValue);
                eFullNamelist.ListBlock[eBlockName].ListBoolValues[eVarName] =
                    eVal;
              } catch (NamelistException &e) {
                parsing_error_end(eBlockName, eVarName, "bool");
              }
            }
            if (eVarNature == "double") {
              double eVal;
              std::istringstream(eVarValue) >> eVal;
              eFullNamelist.ListBlock[eBlockName].ListDoubleValues[eVarName] =
                  eVal;
            }
            if (eVarNature == "string") {
              try {
                std::string eVal =
                    NAMELIST_ConvertFortranStringToCppString(eVarValue);
                eFullNamelist.ListBlock[eBlockName].ListStringValues[eVarName] =
                    eVal;
              } catch (NamelistException &e) {
                parsing_error_end(eBlockName, eVarName, "string");
              }
            }
            if (eVarNature == "listdouble") {
              std::vector<double> eVal =
                  NAMELIST_ConvertFortranStringListDoubleToCppVectorDouble(
                      eVarValue);
              eFullNamelist.ListBlock[eBlockName]
                  .ListListDoubleValues[eVarName] = eVal;
            }
            if (eVarNature == "listint") {
              std::vector<int> eVal =
                  NAMELIST_ConvertFortranStringListIntToCppVectorInt(eVarValue);
              eFullNamelist.ListBlock[eBlockName].ListListIntValues[eVarName] =
                  eVal;
            }
            if (eVarNature == "liststring") {
              try {
                std::vector<std::string> eVal =
                    NAMELIST_ConvertFortranListStringToCppListString(eVarValue);
                eFullNamelist.ListBlock[eBlockName]
                    .ListListStringValues[eVarName] = eVal;
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
}

void NAMELIST_ReadListString(FullNamelist &eFullNamelist, std::vector<std::string> const &ListString) {
  size_t lenString = 30;
  std::string eRandString = random_string(lenString);
  std::string ePrefix = "/tmp/Std_adm";
  std::string TheFile = ePrefix + eRandString;
  {
    std::ofstream OUTfs(TheFile);
    for (auto const &eStr : ListString)
      OUTfs << eStr << "\n";
  }
  NAMELIST_ReadNamelistFile(TheFile, eFullNamelist);
}



// clang-format off
#endif  // SRC_BASIC_NAMELIST_H_
// clang-format on
