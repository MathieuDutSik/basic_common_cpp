// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_file.h"
#include <chrono>
#include <filesystem>

void check_external(bool cond, std::string const &msg) {
  if (!cond) {
    std::cerr << "External program test failure: " << msg << "\n";
    throw TerminalException{1};
  }
}

int main() {
  try {
    namespace fs = std::filesystem;

    auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    fs::path root = fs::temp_directory_path() /
                    ("basic_common_cpp_external_test_" + std::to_string(now));
    fs::remove_all(root);
    fs::create_directories(root);

    std::string out_file = (root / "stdout.txt").string();
    std::string err_file = (root / "stderr.txt").string();
    std::string in_file = (root / "stdin.txt").string();
    std::string cat_out_file = (root / "cat_stdout.txt").string();

    {
      std::ofstream os(in_file);
      os << "first line\n";
      os << "second line\n";
    }

    int code = RunExternalProgram("/bin/echo", {"hello world"}, out_file, std::nullopt);
    check_external(code == 0, "echo should exit with code 0");
    std::vector<std::string> out_lines = ReadFullFile(out_file);
    check_external(out_lines.size() == 1, "echo output should have one line");
    check_external(out_lines[0] == "hello world", "echo output should match argument");

    std::string missing_path = (root / "missing-entry").string();
    code = RunExternalProgram("/bin/ls", {missing_path}, std::nullopt, err_file);
    check_external(code != 0, "ls on a missing path should fail");
    check_external(IsExistingFile(err_file), "stderr redirection file should exist");
    check_external(FILE_GetNumberLine(err_file) > 0, "stderr redirection file should be non-empty");

    code = RunExternalProgramWithInputFile("/bin/cat", {}, in_file, cat_out_file,
                                           std::nullopt);
    check_external(code == 0, "cat with redirected stdin should exit with code 0");
    check_external(ReadFullFile(cat_out_file) == ReadFullFile(in_file),
                   "cat output should match redirected stdin");

    fs::remove_all(root);
    std::cerr << "Normal termination of Test_external_program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_external_program\n";
    exit(e.eVal);
  }
}
