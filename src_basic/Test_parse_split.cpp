// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_string.h"

int main(int argc, char *argv[]) {
  std::string eLine(
      "Reading existing matrix=37 8 4 1 4 1 1 4 1 1 2 4 -2 0 -2 -2 4 0 2 2 2 "
      "-1 4 2 0 1 1 -2 -1 4 2 0 1 1 -2 -1 1 4 idxMatrixCurrent=10893END");
  std::optional<std::vector<std::string>> opt = STRING_ParseSingleLine(
      eLine, {"Reading existing matrix=", " idxMatrixCurrent=", "END"});
  if (!opt) {
    std::cerr << "Failed to parse eLine\n";
    throw TerminalException{1};
  }
  std::vector<std::string> const &LStrParse = *opt;

  std::cerr << "|LStrParse|=" << LStrParse.size() << "\n";
  std::cerr << "LStrParse[0]=" << LStrParse[0] << "\n";
  std::cerr << "LStrParse[1]=" << LStrParse[1] << "\n";
}
