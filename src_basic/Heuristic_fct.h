// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_HEURISTIC_FCT_H_
#define SRC_BASIC_HEURISTIC_FCT_H_

#include "Basic_file.h"
#include "Temp_common.h"
#include "Namelist.h"
#include <map>
#include <string>
#include <vector>

template <typename T> struct SingleCondition {
  std::string eCond;
  std::string eType;
  T NumValue;
};

void CheckEType(std::string const& eType)
{
  std::vector<std::string> LTypes{">", ">=", "=", "<", "<="};
  if (PositionVect(LTypes, eType) == -1) {
    std::cerr << "We found eType=" << eType << "\n";
    std::cerr << "But the allowed types are";
    for (auto & eStr : LTypes)
      std::cerr << " " << eStr;
    std::cerr << "\n";
    throw TerminalException{1};
  }
}


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
      CheckEType(eType);
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
      auto get_value=[&]() -> T {
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
      };
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

template <typename T>
TheHeuristic<T>
HeuristicFrom_LS_LS_S(std::string const &Default,
                      std::vector<std::string> const &ListFullCond,
                      std::vector<std::string> const &ListConclusion) {
  std::vector<OneFullCondition<T>> AllTests;
  if (ListCond.size() != ListConclusion.size()) {
    std::cerr << "ListCond length different from ListConclusion length\n";
    throw TerminalException{1};
  }
  for (size_t i=0; i<ListCond.size(); i++) {
    std::string const& eFullCond = ListFullCond[i];
    size_t len = eFullCond.size();
    std::string char1 = "(";
    std::string char2 = ")";
    for (size_t iC=0; iC<len; iC++) {
      std::string eC = eFullCond.substr(iC,1);
      if (eC == char1 || eC == char2) {
        std::cerr << "Use of forbidden character eC=" << eC << " in the condition\n";
        throw TerminalException{1};
      }
    }
    std::vector<SingleCondition<T>> TheConditions;
    std::vector<std::string> LStr = STRING_Split(eFullCond, "&&");
    for (auto & eStr : LStr) {
      std::vector<std::string> LStrB = STRING_Split(eStr, " ");
      if (LStrB.size() != 3) {
        std::cerr << "|LStrB|=" << LStrB.size() << "\n";
        throw TerminalException{1};
      }
      std::string eCond = LStrB[0];
      std::string eType = LStrB[1];
      CheckEType(eType);
      T eNum = ParseScalar<T>(LStrB[2]);
      SingleCondition<T> eSingCond{eCond, eType, eNum};
      TheConditions.push_back(eSingCond);
    }
    std::string TheResult = ListConclusion[i];
    AllTests.push_back({TheConditions, TheResult});
  }
  return {std::move(AllTests), Default};
}

template <typename T>
void CheckHeuristicInput(TheHeuristic<T> const& heu, std::vector<std::string> const& ListAllowedInput)
{
  for (auto & eFullCond : heu.AllTests) {
    for (auto & eCond : eFullCond.TheConditions) {
      if (PositionVect(ListAllowedInput, eCond) == -1) {
        std::cerr << "The variable eCond=" << eCond << " is not allowed on input\n";
        std::cerr << "ListAlowedInput =";
        for (auto & eStr : ListAllowedInput)
          std::cerr << " " << eStr;
        std::cerr << "\n";
        throw TerminalException{1};
      }
    }
  }
}

template <typename T>
void CheckHeuristicOutput(TheHeuristic<T> const& heu, std::vector<std::string> const& ListAllowedOutput)
{
  auto check=[&](std::string const& output) -> void {
    if (PositionVect(ListAllowedOutput, output) == -1) {
      std::cerr << "The possible output=" << eCond << " is not allowed\n";
      std::cerr << "ListAlowedOutput =";
      for (auto & eStr : ListAllowedOutput)
        std::cerr << " " << eStr;
      std::cerr << "\n";
      throw TerminalException{1};
    }
  };
  for (auto & eFullCond : heu.AllTests) {
    check(eFullCond.TheResult);
  }
  check(eFullCond.DefaultResult);
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


template <typename T>
std::optional<TimingComputationResult<T>> ReadTimingComputationResult(std::string const& line, std::string const& name)
{
  std::optional<std::vector<std::string>> opt;
  opt = STRING_ParseSingleLine(line, {"name=", " keys=(", ") choice=", " result=", " END"});
  if (!opt)
    return {};
  std::vector<std::string> LStr = *opt;
  if (LStr[0] != name)
    return {};
  //
  std::string str_keys = LStr[1];
  std::vector<std::string> LStrKey = STRING_Split(str_keys, ", ");
  std::vector<std::pair<std::string,T>> keys;
  for (auto& eStr : LStr) {
    std::vector<std::string> LStrB = STRING_Split(eStr, ":");
    if (LStrB.size() != 2)
      return {};
    std::string key = LStrB[0];
    std::string str_val = LStrB[1];
    T val = ParseScalar<T>(str_val);
    std::pair<std::string,T> ep{key,val};
    keys.push_back(ep);
  }
  //
  std::string choice = LStr[2];
  //
  std::string str_result = LStr[3];
  size_t result = ParseScalar<size_t>(str_result);
  //
  TimingComputationAttempt<T> tca{keys, choice};
  TimingComputationResult<T> tcr{tca,result};
  return tcr;
}

template <typename T>
void PrintTimingComputationResult(std::ostream& os, TimingComputationResult<T> const& eTCR, std::string const& name)
{
  os << "name=" << name;
  os << " keys=(";
  bool IsFirst = true;
  for (auto & eEnt : eTCR.input.keys) {
    if (!IsFirst) {
      os << ", ";
    }
    IsFirst = false;
    os << eEnt.first << ":" << eEnt.second;
  }
  os << ") ";
  os << " choice=" << eTCR.input.choice;
  os << " result=" << eTCR.result;
  os << " END\n";
}


struct LimitedEmpiricalDistributionFunction {
  size_t n_max;
  size_t n_ins;
  std::map<double,size_t> ListValWei;
  LimitedEmpiricalDistributionFunction(size_t n_max) : n_max(n_max), n_ins(0) {
  }
  void clear_entry() {
    //    std::cerr << "|ListValWei|=" << ListValWei.size() << "\n";
    auto iter = ListValWei.begin();
    double min_delta = std::numeric_limits<double>::max();
    size_t pos = 0;
    size_t pos_found = 0;
    // Determine the smallest delta
    while(true) {
      auto iterN = iter;
      iterN++;
      if (iterN == ListValWei.end())
        break;
      //      std::cerr << "pos=" << pos << "  val=" << iter->first << " mult=" << iter->second << "\n";
      double delta = iterN->first - iter->first;
      //      std::cerr << "   delta=" << delta << "\n";
      if (delta < min_delta) {
        pos_found = pos;
        min_delta = delta;
      }
      pos++;
      iter++;
    }
    //    std::cerr << "NEXT : pos_found=" << pos_found << " |ListValWei|=" << ListValWei.size() << "\n";
    // Now clearing the entries
    auto iter1 = ListValWei.begin();
    for (size_t u=0; u<pos_found; u++)
      iter1++;
    auto iter2 = iter1;
    iter2++;
    double val1 = iter1->first;
    double val2 = iter2->first;
    double w1 = iter1->second;
    double w2 = iter2->second;
    double new_val = (val1 * w1 + val2 * w2) / (w1 + w2);
    double new_w = w1 + w2;
    ListValWei.erase(val1);
    ListValWei.erase(val2);
    ListValWei[new_val] = new_w;
  }
  void insert_value(double new_val) {
    ListValWei[new_val] += 1;
    n_ins++;
    if (ListValWei.size() > n_max) {
      clear_entry();
    }
  }
  double get_percentile(double const& alpha) const {
    size_t sum_w = 0;
    auto iter = ListValWei.begin();
    size_t crit_w = round(alpha * n_ins);
    while(true) {
      sum_w += iter->second;
      if (sum_w > crit_w)
        return iter->first;
      iter++;
      if (iter == ListValWei.end()) {
        std::cerr << "Failed to find an entry in the map\n";
        throw TerminalException{1};
      }
    }
  }
  std::string string() const {
    std::string str_ret;
    bool IsFirst = true;
    auto iter = ListValWei.begin();
    while(true) {
      if (!IsFirst)
        str_ret += ", ";
      if (IsFirst)
        IsFirst = false;
      str_ret += "(" + std::to_string(iter->first) + "," + std::to_string(iter->second) + ")";
      iter++;
      if (iter == ListValWei.end()) {
        break;
      }
    }
    return str_ret;
  }
};


FullNamelist NAMELIST_GetStandard_RecursiveDualDescription() {
  std::map<std::string, SingleBlock> ListBlock;
  // PROBABILITY DISTRIBUTIONS
  {
    std::map<std::string, std::vector<std::string>> ListListStringValues;
    std::map<std::string, std::vector<int>> ListListIntValues;
    ListListStringValues["ListName"] = {"distri1", "distri2", "distri3", "distri4"};
    ListListIntValues["ListN"] = {100, 100, -1, -1};
    ListListStringValues["ListNature"] = {"dirac", "sampled", "file_select", "file"};
    ListListStringValues["ListDescription"] = {"145", "130:50 145:30 150:20", "InputA_key1_val1_key2_val2", "InputB"};
    SingleBlock BlockPROBA;
    BlockPROBA.ListListStringValues = ListListStringValues;
    BlockPROBA.ListListIntValues = ListListIntValues;
    ListBlock["PROBABILITY_DISTRIBUTIONS"] = BlockPROBA;
  }
  // SINGLE THOMPSON STATE
  {
    std::map<std::string, std::vector<std::string>> ListListStringValues;
    ListListStringValues["ListAnswer"] = {"cdd", "lrs"};
    ListListStringValues["ListName"] = {"state1", "state2"};
    ListListStringValues["ListDescription"] = {"cdd:distri1 lrs:distri2", "cdd:distri1 lrs:distri3"};
    SingleBlock BlockTHOMPSON;
    BlockTHOMPSON.ListListStringValues = ListListStringValues;
    ListBlock["THOMPSON_STATE"] = BlockTHOMPSON;
  }
  // GROUPING
  {
    std::map<std::string, std::vector<std::string>> ListListStringValues;
    ListListStringValues["ListKey"] = {"incidence", "groupsize"};
    ListListStringValues["ListDescription"] = {"superfine", "1-10,11-20,21-infinit"};
    SingleBlock BlockGROUP;
    BlockGROUP.ListListStringValues = ListListStringValues;
    ListBlock["GROUPINGS"] = BlockGROUP;
  }
  // PRIOR
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, bool> ListBoolValues;
    std::map<std::string, std::vector<std::string>> ListListStringValues;
    ListStringValues["DefaultPrior"] = "noprior";
    ListListStringValues["ListFullCond"] = {"incidence > 15 && groupsize > 10", "incidence GroupEXT"};
    ListListStringValues["ListConclusion"] = {"state1", "state2"};
    SingleBlock BlockPRIOR;
    BlockPRIOR.ListBoolValues = ListBoolValues;
    BlockPRIOR.ListStringValues = ListStringValues;
    BlockPRIOR.ListListStringValues = ListListStringValues;
    ListBlock["PRIORS"] = BlockPRIOR;
  }
  // IO
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, bool> ListBoolValues;
    ListStringValues["name"] = "split";
    ListBoolValues["ProcessExistingDataIfExist"] = false;
    ListBoolValues["WriteLog"] = false;
    ListStringValues["LogFileProcessing"] = "irrelevant";
    SingleBlock BlockIO;
    BlockIO.ListBoolValues = ListBoolValues;
    BlockIO.ListStringValues = ListStringValues;
    ListBlock["IO"] = BlockIO;
  }
  // Merging all data
  return {std::move(ListBlock), "undefined"};
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
//
// Now the algorithm itself.
// --- We need to keep the possibility of bypassing the Thompson sampling
//   when we so want.
//   This can be done by having an attribute FORC to some of the rules.
// --- From the initial heuristic we can get many information:
//    --- Which variables to consider.
//    --- The relative importance of the variable. First the number of
//     vertices, then the size of the group.
//    --- What are the allowed outputs.
//    --- The maximal effective range of variables. It is true that for
//      things like group size, an exponential profile may be better.
//      but for the decision, that should actually be neutral. Since group
//      size that triggers are typically small (say 100) in contrast to smaller
// --- We can store the measurements and from that do some interpolation.
// --- We can decrease the number of measurements used if it gets to wild.
//   This can be done by merging points.
// --- For interpolation, taking the inverse sqare distance is probably the
//   best.
//
// The Thompson sampling goes that way:
// ---We have a number of measurements (or probability distribution) for each
//    methods and inputs.
// ---For each we compute the expected outcome of a choice by sampling the
//    probability distribution.
// ---We insert the result into the list of measurement (or update the probability
//    distribution).
// If we are using probability distributions, then it should have nice probability
// distributions, see https://en.wikipedia.org/wiki/Conjugate_prior
//
// Further thinking:
// ---We can think of the initial heuristic as a kind of prior put for the heuristic.
// ---This solves following problem:
//    ---If we want to force the decision in some case, we simply put the distribution
//       appropriately.
//    ---It allows to integrate the heuristics with the Thompson sampling.
//    ---
// ---This cause many major problems:
//    ---It forces us to put numbers for the runtime.
//    ---We have to set up some probability distributions for each case and that
//       can be intensive in term of text files and thinking about the input.
//    ---We need to choose the sampling size at the beginning.
//    ---We need to be able to start on no-prior, all methods equally acceptable at the start.
//    ---We can write a
//    ---We can two many categories
// ---The problem is that this force us to give runtimes in advance. And this can
//    be an issue since we have to effectively measure by ourselves
// ---We can have heuristic grouped by 
//
template <typename T>
struct ThompsonSamplingHeuristic {
  std::ostream& os;
  std::string name;
  bool WriteLog;
  TheHeuristic<T> heu;
  std::vector<TimingComputationResult<T>> l_completed_info;
  std::vector<TimingComputationAttempt<T>> l_submission;
  ThompsonSamplingHeuristic(std::ostream& os, FullNamelist const& TS) : os(os) {
    SingleBlock const& BlockPROBA = TS.ListBlock.at("PROBABILITY_DISTRIBUTIONS");
    SingleBlock const& BlockTHOMPSON = TS.ListBlock.at("THOMPSON_STATE");
    SingleBlock const& BlockPRIOR = TS.ListBlock.at("PRIORS");
    SingleBlock const& BlockIO = TS.ListBlock.at("IO");
    name = BlockIO.ListStringValues.at("name");
    WriteLog = BlockIO.ListBoolValues.at("WriteLog");
    //
    std::string const& DefaultPrior = BlockPRIOR.ListStringValues.at("DefaultPrior");
    std::vector<std::string> const& ListFullCond = BlockPRIOR.ListListStringValues.at("ListFullCond");
    std::vector<std::string> const& ListConclusion = BlockPRIOR.ListListStringValues.at("ListConclusion");
    heu = HeuristicFrom_LS_LS_S<T>(DefaultPrior, ListFullCond, ListConclusion);

  }
  void InsertCompletedInfo(std::string const& file) {
    if (!IsExistingFile(file)) {
      std::cerr << "The file=" << file << " is missing.\n";
      std::cerr << "Cannot parse it in SelfCorrectingHeuristic with name=" << name << "\n";
      throw TerminalException{1};
    }
    std::ifstream is(file);
    std::string line;
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
    TimingComputationAttempt<T> einput = l_submission.front();
    TimingComputationResult<T> eTCR{einput, result};
    l_completed_info.push_back(eTCR);
    l_submission.pop();
    if (WriteLog) {
      PrintTimingComputationResult(os, eTCR, name);
    }
  }
  ~ThompsonSamplingHeuristic() {
    if (l_submission.size() > 0) {
      std::cerr << "The SelfCorrectingHeuristic terminates with some submission being uncompleted\n";
      std::cerr << "Just so you know. Most likely, premature termination, otherwise a bug\n";
    }
  }
};





// clang-format off
#endif  // SRC_BASIC_HEURISTIC_FCT_H_
// clang-format on
