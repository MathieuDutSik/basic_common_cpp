// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_EXTERNAL_PROGRAM_WINDOWS_H_
#define SRC_BASIC_BASIC_EXTERNAL_PROGRAM_WINDOWS_H_

// Windows / MinGW implementation: there is no fork / execvp / waitpid, so we use
// the CRT process facilities instead. _spawnvp(_P_WAIT, ...) runs the program
// synchronously and returns its exit code, and the std_input / std_output /
// std_error redirection is done by temporarily replacing fds 0 / 1 / 2 with the
// requested files (the spawned child inherits them) and restoring them
// afterwards. Included through Basic_external_program.h; do not include this
// header directly.

#include "Temp_common.h"
#include <fcntl.h>
#include <io.h>
#include <optional>
#include <process.h>
#include <string>
#include <sys/stat.h>
#include <vector>

int RunExternalProgram(std::string const &program,
                       std::vector<std::string> const &inputs,
                       std::optional<std::string> const &std_input,
                       std::optional<std::string> const &std_output,
                       std::optional<std::string> const &std_error) {
  if (program.empty()) {
    std::cerr << "RunExternalProgram requires a non-empty program name\n";
    throw TerminalException{1};
  }
  int saved_in = -1;
  int saved_out = -1;
  int saved_err = -1;
  auto restore_fds = [&]() {
    if (saved_in != -1) {
      _dup2(saved_in, 0);
      _close(saved_in);
      saved_in = -1;
    }
    if (saved_out != -1) {
      _dup2(saved_out, 1);
      _close(saved_out);
      saved_out = -1;
    }
    if (saved_err != -1) {
      _dup2(saved_err, 2);
      _close(saved_err);
      saved_err = -1;
    }
  };
  if (std_input) {
    int fd = _open(std_input->c_str(), _O_RDONLY | _O_BINARY);
    if (fd == -1) {
      std::cerr << "RunExternalProgram failed to open input file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "std_input=" << *std_input << "\n";
      throw TerminalException{1};
    }
    saved_in = _dup(0);
    _dup2(fd, 0);
    _close(fd);
  }
  if (std_output) {
    int fd = _open(std_output->c_str(),
                   _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
                   _S_IREAD | _S_IWRITE);
    if (fd == -1) {
      restore_fds();
      std::cerr << "RunExternalProgram failed to open output file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "std_output=" << *std_output << "\n";
      throw TerminalException{1};
    }
    saved_out = _dup(1);
    _dup2(fd, 1);
    _close(fd);
  }
  if (std_error) {
    int fd = _open(std_error->c_str(),
                   _O_WRONLY | _O_CREAT | _O_TRUNC | _O_BINARY,
                   _S_IREAD | _S_IWRITE);
    if (fd == -1) {
      restore_fds();
      std::cerr << "RunExternalProgram failed to open error file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "std_error=" << *std_error << "\n";
      throw TerminalException{1};
    }
    saved_err = _dup(2);
    _dup2(fd, 2);
    _close(fd);
  }
  std::vector<std::string> args;
  args.reserve(inputs.size() + 1);
  args.push_back(program);
  for (auto const &arg : inputs) {
    args.push_back(arg);
  }
  std::vector<char *> argv;
  argv.reserve(args.size() + 1);
  for (auto &arg : args) {
    argv.push_back(arg.data());
  }
  argv.push_back(nullptr);
  intptr_t rc = _spawnvp(_P_WAIT, program.c_str(), argv.data());
  restore_fds();
  if (rc == -1) {
    std::cerr << "RunExternalProgram failed in _spawnvp\n";
    std::cerr << "program=" << program << "\n";
    throw TerminalException{1};
  }
  return static_cast<int>(rc);
}

// clang-format off
#endif  // SRC_BASIC_BASIC_EXTERNAL_PROGRAM_WINDOWS_H_
// clang-format on
