// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_HEURISTIC_FCT_H_
#define SRC_BASIC_HEURISTIC_FCT_H_

#include "Basic_file.h"
#include "Temp_common.h"
#include <map>
#include <string>
#include <vector>

template <typename T> struct SingleCondition {
  std::string eCond;
  std::string eType;
  T NumValue;
};

template <typename T>
std::ostream &operator<<(std::ostream &os,
                         SingleCondition<T> const &eSingCond) {
  os << eSingCond.eCond << " " << eSingCond.eType << " " << eSingCond.NumValue;
  return os;
}

template <typename T> struct OneFullCondition {
  std::vector<SingleCondition<T>> TheConditions;
  std::string TheResult;
};

template <typename T>
std::ostream &operator<<(std::ostream &os,
                         OneFullCondition<T> const &eFullCond) {
  size_t len = eFullCond.TheConditions.size();
  for (size_t i = 0; i < len; i++) {
    if (i > 0)
      os << " && ";
    os << "(" << eFullCond.TheConditions[i] << ")";
  }
  os << " => " << eFullCond.TheResult;
  return os;
}

template <typename T> struct TheHeuristic {
  std::vector<OneFullCondition<T>> AllTests;
  std::string DefaultResult;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, TheHeuristic<T> const &Heu) {
  size_t len = Heu.AllTests.size();
  for (size_t i = 0; i < Heu.AllTests.size(); i++) {
    os << "   i=" << i << "/" << len << " " << Heu.AllTests[i] << "\n";
  }
  os << "      Default=" << Heu.DefaultResult;
  return os;
}

template <typename T> TheHeuristic<T> ReadHeuristic(std::istream &is) {
  if (!is.good()) {
    std::cerr
        << "ReadHeuristic operation failed because stream is not valied\n";
    throw TerminalException{1};
  }
  TheHeuristic<T> TheHeu;
  int nbFullCond;
  is >> nbFullCond;
  if (nbFullCond < 0) {
    std::cerr << "We must have nbFullCond >= 0\n";
    throw TerminalException{1};
  }
  for (int iFullCond = 0; iFullCond < nbFullCond; iFullCond++) {
    std::vector<SingleCondition<T>> ListSingleCond;
    int nbCond;
    is >> nbCond;
    if (nbCond <= 0) {
      std::cerr << "iFullCond=" << iFullCond << " nbFullCond=" << nbFullCond
                << "\n";
      std::cerr << "Error, we must have nbCond > 0\n";
      std::cerr << "nbCond=" << nbCond << "\n";
      throw TerminalException{1};
    }
    for (int iCond = 0; iCond < nbCond; iCond++) {
      std::string eType, eCond;
      T eNum;
      is >> eCond;
      is >> eType;
      is >> eNum;
      std::vector<std::string> ListType{">", "<", "=", "<=", ">="};
      bool IsMatch = false;
      for (auto &eTypePos : ListType) {
        if (eTypePos == eType)
          IsMatch = true;
      }
      if (!IsMatch) {
        std::cerr
            << "Only allowed possibilities for eType are <, >, =, <= and >=\n";
        throw TerminalException{1};
      }
      if (eCond.size() == 0) {
        std::cerr << "eCond must be nontrivial otherwise evaluation will be "
                     "impossible\n";
        throw TerminalException{1};
      }
      SingleCondition<T> eSingCond{eCond, eType, eNum};
      ListSingleCond.push_back(eSingCond);
    }
    std::string eResult;
    is >> eResult;
    if (eResult.size() == 0) {
      std::cerr << "eResult must be nontrivial otherwise evaluation will be "
                   "impossible\n";
      throw TerminalException{1};
    }
    OneFullCondition<T> OneCond{ListSingleCond, eResult};
    TheHeu.AllTests.push_back(OneCond);
  }
  std::string DefaultResult;
  is >> DefaultResult;
  if (DefaultResult.size() == 0) {
    std::cerr << "DefaultResult must be nontrivial otherwise evaluation will "
                 "be impossible\n";
    throw TerminalException{1};
  }
  TheHeu.DefaultResult = DefaultResult;
  return TheHeu;
}

template <typename T>
void ReadHeuristicFileCond(std::string const &eFile, TheHeuristic<T> &eHeu) {
  if (eFile != "unset.heu") {
    std::cerr << "eFile=" << eFile << "\n";
    IsExistingFileDie(eFile);
    std::ifstream is(eFile);
    try {
      eHeu = ReadHeuristic<T>(is);
    } catch (TerminalException const &e) {
      std::cerr << "Failed in reading the file eFile=" << eFile << "\n";
      throw TerminalException{1};
    }
  }
}

template <typename T>
std::string HeuristicEvaluation(std::map<std::string, T> const &TheCand,
                                TheHeuristic<T> const &TheHeu) {
  for (auto const &eFullCond : TheHeu.AllTests) {
    bool IsOK = true;
    for (auto const &eSingCond : eFullCond.TheConditions) {
      std::string eCond = eSingCond.eCond;
      std::string eType = eSingCond.eType;
      T eNum = eSingCond.NumValue;
      auto get_value[&]() -> T {
        auto search = TheCand.find(eCond);
        if (search != TheCand.end())
          return search->second;
        std::cerr << "Entry " << eCond << " is required by heuristic\n";
        std::cerr << "Yet it is missing in the Candidate. TheCand=\n";
        for (auto &kv : TheCand) {
          std::cerr << "  key=" << kv.first << " value=" << kv.second << "\n";
        }
        std::cerr << "Please correct\n";
        throw TerminalException{1};
      }
      T eValue = get_value();
      bool WeMatch = false;
      if (eValue > eNum && eType == ">")
        WeMatch = true;
      if (eValue >= eNum && eType == ">=")
        WeMatch = true;
      if (eValue == eNum && eType == "=")
        WeMatch = true;
      if (eValue < eNum && eType == "<")
        WeMatch = true;
      if (eValue <= eNum && eType == "<=")
        WeMatch = true;
      if (!WeMatch)
        IsOK = false;
    }
    if (IsOK)
      return eFullCond.TheResult;
  }
  return TheHeu.DefaultResult;
}

template <typename T>
TheHeuristic<T>
HeuristicFromListString(std::vector<std::string> const &ListString) {
  size_t lenString = 30;
  std::string eRandString = random_string(lenString);
  std::string ePrefix = "/tmp/Std_adm";
  std::string TheFile = ePrefix + eRandString;
  std::ofstream OUTfs(TheFile);
  for (auto const &eStr : ListString)
    OUTfs << eStr << "\n";
  OUTfs.close();
  // Now reading it
  std::ifstream INfs(TheFile);
  TheHeuristic<T> TheHeu = ReadHeuristic<T>(INfs);
  std::remove(TheFile.c_str());
  return TheHeu;
}


// Standing for TCS = "Timing Computation Result"
template <typename T>
struct TimingComputationAttempt {
  std::vector<std::pair<std::string,T>> keys;
  std::string choice;
};



template <typename T>
struct TimingComputationResult {
  TimingComputationAttempt<T> input;
  size_t result;
};

std::vector<std::pair<std::string,std::string>> ParseLineSeparator(std::string const& line)
{
  
}


std::optional<TimingComputationResult> ReadTimingComputationResult(std::string const& line, std::string const& name)
{
  
}

void PrintTimingComputationResult(std::ostream& os, TimingComputationResult const& eTCS, std::string const& name)
{
  os << "name=" << name << " ";
  os << "keys=(";
  bool IsFirst = true;
  for (auto & eEnt : eTCS.input.keys) {
    if (!IsFirst) {
      os << ", ";
    }
    IsFirst = false;
    os << eEnt.first << ":" << eEnt.second;
  }
  os << ") ";
  os << " choice=" << eTCS.input.choice << " ";
  os << "result=" << eTCS.result;
  os << "\n";
}



//
// We want to improve the heuristic stuff so that it is handled
// directly by the code
//
// Design:
// --- Things are done in two steps:
//    --- First submitting the query and getting the result
//    --- Then getting the result (and printing it)
// --- For queries like "split", we have multiple queries going on
//   where we need to wait for the result to be obtained.
//   Therefore, this forces on us a recursive designs where we accumulate
//   the queries one by one and then update later on.
// --- Would it be better to have one data structure for "split", one for
//   "DualDesc", one for XXX instead of just one bg data structure?
//   Separate data structures are probably the better approach here.
//   This is forced on us for the reading of information from file.
// --- Actually a better design is to have a "name" for the heuristic so
//   that it knows where to read the data.
template <typename T>
struct SelfCorrectingHeuristic {
  std::string name;
  std::ostream& os;
  bool DoPrint;
  TheHeuristic<T> heu;
  std::vector<TimingComputationResult<T>> l_completed_info;
  std::vector<TimingComputationAttempt<T>> l_submission;
  
  
  SelfCorrectingHeuristic(std::string const& name, std::ostream& os, bool DoPrint) : name(name), os(_os), DoPrint(DoPrint) {
  }
  void InsertCompletedInfo(std::string const& file) {
    if (!IsExistingFile(file)) {
      std::cerr << "The file=" << file << " is missing.\n";
      std::cerr << "Cannot parse it in SelfCorrectingHeuristic with name=" << name << "\n";
      throw TerminalException{1};
    }
    std::istream is(file);
    while (std::getline(is, line)) {
      std::optional<TimingComputationResult<T>> opt = ReadTimingComputationResult<T>(line, name);
      if (opt) {
        l_completed_info.push_back(*opt);
      }
    }
  }
  std::string Kernel_GetEvaluation(std::map<std::string,T> const& TheCand) {
    // Need to add the stuff
    return HeuristicEvaluation(TheCand,heu);
  }
  std::string GetEvaluation(std::map<std::string,T> const& TheCand) {
    std::string choice = Kernel_GetEvaluation(TheCand);
    std::vector<std::pair<std::string,T>> keys;
    for (auto & kv : TheCand) {
      std::pair<std::string,T> ep{kv.first, kv.second};
      keys.push_back(ep);
    }
    TimingComputationAttempt<T> tca{keys, choice};
    l_submission.push_back(tca);
  }
  void SubmitResult(size_t result) {
    Tinput einput = l_submission.front();
    TimingComputationResult eTCS{einput, result};
    l_completed_info.push_back(eTCS);
    l_submission.pop();
    if (DoPrint) {
      PrintTimingComputationResult(os, eTCS, name);
    }
  }
  ~SelfCorrectingHeuristic() {
    if (l_submission.size() > 0) {
      std::cerr << "The SelfCorrectingHeuristic terminates with some submission being uncompleted\n";
      std::cerr << "Just so you know. Most likely, premature termination, otherwise a bug\n";
    }
  }
}





// clang-format off
#endif  // SRC_BASIC_HEURISTIC_FCT_H_
// clang-format on
