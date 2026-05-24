// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Heuristic_ThompsonSampling.h"

using T = int;

static void check(bool cond, std::string const &msg) {
  if (!cond) {
    std::cerr << "Test failed: " << msg << "\n";
    throw TerminalException{1};
  }
}

template <typename F>
static void expect_throw(F const &f, std::string const &msg) {
  bool caught = false;
  try {
    f();
  } catch (TerminalException const &) {
    caught = true;
  }
  if (!caught) {
    std::cerr << "Expected throw but none happened: " << msg << "\n";
    throw TerminalException{1};
  }
}

// Two rules over keys {incidence, groupsize}:
//   (incidence > 10 && groupsize >= 5) => cdd
//   (incidence < 5)                    => lrs
//   default                            => glrs
static TheHeuristic<T> make_sample_heuristic() {
  std::vector<std::string> ListString = {
      "2",
      "2 incidence > 10 groupsize >= 5 cdd",
      "1 incidence < 5 lrs",
      "glrs"};
  return HeuristicFromListString<T>(ListString);
}

static std::string eval(TheHeuristic<T> const &heu, T inc, T grp) {
  std::map<std::string, T> cand{{"incidence", inc}, {"groupsize", grp}};
  return HeuristicEvaluation(cand, heu);
}

static void test_evaluation_default_and_matches() {
  TheHeuristic<T> heu = make_sample_heuristic();
  check(eval(heu, 7, 3) == "glrs", "no rule matches -> default");
  check(eval(heu, 11, 5) == "cdd", "first rule matches");
  check(eval(heu, 4, 0) == "lrs", "second rule matches");
}

static void test_evaluation_boundary() {
  TheHeuristic<T> heu = make_sample_heuristic();
  // 10 is NOT > 10, and not < 5, so default
  check(eval(heu, 10, 10) == "glrs", "strict > boundary should not match");
  // groupsize == 5 is >= 5, with incidence > 10 first rule matches
  check(eval(heu, 100, 5) == "cdd", "non-strict >= boundary should match");
  // first rule wins even if second would also conceptually apply: build a case
  // where both rules cannot both fire (they are mutually exclusive on
  // incidence), so confirm first-match-wins by adding overlap below.
}

static void test_first_match_wins() {
  // Rule A: (x > 0) => A, Rule B: (x > 5) => B
  // For x=10 both match; A is listed first and must win.
  std::vector<std::string> ListString = {
      "2",
      "1 x > 0 A",
      "1 x > 5 B",
      "D"};
  TheHeuristic<T> heu = HeuristicFromListString<T>(ListString);
  std::map<std::string, T> cand{{"x", 10}};
  check(HeuristicEvaluation(cand, heu) == "A", "first matching rule wins");
}

static void test_each_operator() {
  auto build = [](std::string const &op, T num,
                  std::string const &result) -> TheHeuristic<T> {
    std::vector<std::string> LStr = {
        "1", "1 x " + op + " " + std::to_string(num) + " " + result,
        "default"};
    return HeuristicFromListString<T>(LStr);
  };
  auto run = [](TheHeuristic<T> const &heu, T v) -> std::string {
    std::map<std::string, T> cand{{"x", v}};
    return HeuristicEvaluation(cand, heu);
  };
  // =
  TheHeuristic<T> heu_eq = build("=", 5, "match");
  check(run(heu_eq, 5) == "match", "x = 5, v=5");
  check(run(heu_eq, 4) == "default", "x = 5, v=4");
  check(run(heu_eq, 6) == "default", "x = 5, v=6");
  // >
  TheHeuristic<T> heu_gt = build(">", 5, "match");
  check(run(heu_gt, 6) == "match", "x > 5, v=6");
  check(run(heu_gt, 5) == "default", "x > 5, v=5");
  check(run(heu_gt, 4) == "default", "x > 5, v=4");
  // >=
  TheHeuristic<T> heu_gte = build(">=", 5, "match");
  check(run(heu_gte, 6) == "match", "x >= 5, v=6");
  check(run(heu_gte, 5) == "match", "x >= 5, v=5");
  check(run(heu_gte, 4) == "default", "x >= 5, v=4");
  // <
  TheHeuristic<T> heu_lt = build("<", 5, "match");
  check(run(heu_lt, 4) == "match", "x < 5, v=4");
  check(run(heu_lt, 5) == "default", "x < 5, v=5");
  check(run(heu_lt, 6) == "default", "x < 5, v=6");
  // <=
  TheHeuristic<T> heu_lte = build("<=", 5, "match");
  check(run(heu_lte, 4) == "match", "x <= 5, v=4");
  check(run(heu_lte, 5) == "match", "x <= 5, v=5");
  check(run(heu_lte, 6) == "default", "x <= 5, v=6");
}

static void test_missing_key_in_candidate() {
  TheHeuristic<T> heu = make_sample_heuristic();
  // Heuristic references "incidence" but candidate has only "groupsize"
  std::map<std::string, T> cand{{"groupsize", 5}};
  expect_throw([&]() { HeuristicEvaluation(cand, heu); },
               "missing key in candidate should throw");
}

static void test_check_etype() {
  CheckEType(">");
  CheckEType(">=");
  CheckEType("=");
  CheckEType("<");
  CheckEType("<=");
  expect_throw([]() { CheckEType("!="); }, "CheckEType !=");
  expect_throw([]() { CheckEType("=="); }, "CheckEType ==");
  expect_throw([]() { CheckEType(""); }, "CheckEType empty");
  expect_throw([]() { CheckEType("foo"); }, "CheckEType foo");
}

static void test_get_input_output() {
  TheHeuristic<T> heu = make_sample_heuristic();
  std::vector<std::string> inputs = GetHeuristicInput(heu);
  // std::set ordering: groupsize < incidence
  check(inputs.size() == 2, "input count");
  check(inputs[0] == "groupsize" && inputs[1] == "incidence",
        "input names sorted");
  std::vector<std::string> outputs = GetHeuristicOutput(heu);
  // Includes the default result too: cdd, glrs, lrs
  check(outputs.size() == 3, "output count");
  check(outputs[0] == "cdd" && outputs[1] == "glrs" && outputs[2] == "lrs",
        "output names sorted, default included");
}

static void test_get_pivots() {
  TheHeuristic<T> heu = make_sample_heuristic();
  std::vector<T> pivots = GetHeuristicPivots(heu, "incidence");
  // Pivots for incidence are derived from values 10 and 5, each producing
  // {v-1, v, v+1}: union = {4, 5, 6, 9, 10, 11}
  std::vector<T> expected{4, 5, 6, 9, 10, 11};
  check(pivots == expected, "pivots for incidence");
  // Key that does not appear -> empty pivot list
  std::vector<T> none = GetHeuristicPivots(heu, "absent");
  check(none.empty(), "pivots for absent key is empty");
}

static void test_check_input_output() {
  TheHeuristic<T> heu = make_sample_heuristic();
  // Allowed set may be a superset
  CheckHeuristicInput(heu, {"incidence", "groupsize", "extra"});
  expect_throw([&]() { CheckHeuristicInput(heu, {"incidence"}); },
               "missing input groupsize should throw");
  CheckHeuristicOutput(heu, {"cdd", "lrs", "glrs", "extra"});
  expect_throw([&]() { CheckHeuristicOutput(heu, {"cdd", "lrs"}); },
               "missing default in allowed outputs should throw");
  expect_throw([&]() { CheckHeuristicOutput(heu, {"glrs", "lrs"}); },
               "missing rule result in allowed outputs should throw");
}

static void test_read_invalid() {
  {
    std::istringstream is("-1\ndefault\n");
    expect_throw([&]() { ReadHeuristic<T>(is); }, "negative nbFullCond");
  }
  {
    std::istringstream is("1\n0\ndefault\n");
    expect_throw([&]() { ReadHeuristic<T>(is); }, "nbCond <= 0");
  }
  {
    std::istringstream is("1\n1\nincidence != 5\ncdd\ndefault\n");
    expect_throw([&]() { ReadHeuristic<T>(is); }, "invalid operator");
  }
}

static void test_read_roundtrip() {
  // The textual stream produced from the same data as make_sample_heuristic
  // must give back an equivalent heuristic.
  std::string text =
      "2\n"
      "2 incidence > 10 groupsize >= 5 cdd\n"
      "1 incidence < 5 lrs\n"
      "glrs\n";
  std::istringstream is(text);
  TheHeuristic<T> heu = ReadHeuristic<T>(is);
  check(eval(heu, 11, 5) == "cdd", "read roundtrip: first rule");
  check(eval(heu, 4, 0) == "lrs", "read roundtrip: second rule");
  check(eval(heu, 7, 0) == "glrs", "read roundtrip: default");
}

static void test_heuristic_from_LS_LS_S() {
  std::vector<std::string> conds = {"incidence > 10 && groupsize >= 5",
                                    "incidence < 5"};
  std::vector<std::string> results = {"cdd", "lrs"};
  TheHeuristic<T> heu = HeuristicFrom_LS_LS_S<T>("glrs", conds, results);
  check(eval(heu, 11, 5) == "cdd", "LS_LS_S: first rule");
  check(eval(heu, 4, 0) == "lrs", "LS_LS_S: second rule");
  check(eval(heu, 7, 0) == "glrs", "LS_LS_S: default");
  // Mismatched conds/results lengths must throw
  expect_throw(
      [&]() {
        HeuristicFrom_LS_LS_S<T>("d", {"x > 0", "y < 1"}, {"only_one"});
      },
      "LS_LS_S: mismatched lengths");
  // Forbidden parenthesis character must throw
  expect_throw(
      [&]() {
        HeuristicFrom_LS_LS_S<T>("d", {"(x > 0)"}, {"A"});
      },
      "LS_LS_S: forbidden character");
  // Malformed single condition (not 3 tokens) must throw
  expect_throw(
      [&]() {
        HeuristicFrom_LS_LS_S<T>("d", {"x >"}, {"A"});
      },
      "LS_LS_S: malformed condition");
}

int main() {
  HumanTime time;
  try {
    test_evaluation_default_and_matches();
    test_evaluation_boundary();
    test_first_match_wins();
    test_each_operator();
    test_missing_key_in_candidate();
    test_check_etype();
    test_get_input_output();
    test_get_pivots();
    test_check_input_output();
    test_read_invalid();
    test_read_roundtrip();
    test_heuristic_from_LS_LS_S();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
