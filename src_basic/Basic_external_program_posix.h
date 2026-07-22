// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_EXTERNAL_PROGRAM_POSIX_H_
#define SRC_BASIC_BASIC_EXTERNAL_PROGRAM_POSIX_H_

// POSIX implementation: runs an external program in a child process via
// fork / execvp / waitpid. Included through Basic_external_program.h; do not
// include this header directly.

#include "Temp_common.h"
#include <fcntl.h>
#include <optional>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
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
  int fd_in = -1;
  int fd_out = -1;
  int fd_err = -1;
  if (std_input) {
    fd_in = open(std_input->c_str(), O_RDONLY);
    if (fd_in == -1) {
      std::cerr << "RunExternalProgram failed to open input file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "std_input=" << *std_input << "\n";
      throw TerminalException{1};
    }
  }
  if (std_output) {
    fd_out = open(std_output->c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_out == -1) {
      if (fd_in != -1) {
        close(fd_in);
      }
      std::cerr << "RunExternalProgram failed to open output file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "std_output=" << *std_output << "\n";
      throw TerminalException{1};
    }
  }
  if (std_error) {
    fd_err = open(std_error->c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_err == -1) {
      if (fd_in != -1) {
        close(fd_in);
      }
      if (fd_out != -1) {
        close(fd_out);
      }
      std::cerr << "RunExternalProgram failed to open error file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "std_error=" << *std_error << "\n";
      throw TerminalException{1};
    }
  }
  pid_t pid = fork();
  if (pid == -1) {
    if (fd_in != -1) {
      close(fd_in);
    }
    if (fd_out != -1) {
      close(fd_out);
    }
    if (fd_err != -1) {
      close(fd_err);
    }
    std::cerr << "RunExternalProgram failed in fork\n";
    std::cerr << "program=" << program << "\n";
    throw TerminalException{1};
  }
  if (pid == 0) {
    if (fd_in != -1) {
      if (dup2(fd_in, STDIN_FILENO) == -1) {
        _exit(127);
      }
    }
    if (fd_out != -1) {
      if (dup2(fd_out, STDOUT_FILENO) == -1) {
        _exit(127);
      }
    }
    if (fd_err != -1) {
      if (dup2(fd_err, STDERR_FILENO) == -1) {
        _exit(127);
      }
    }
    if (fd_in != -1) {
      close(fd_in);
    }
    if (fd_out != -1) {
      close(fd_out);
    }
    if (fd_err != -1) {
      close(fd_err);
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
    execvp(program.c_str(), argv.data());
    _exit(127);
  }
  if (fd_in != -1) {
    close(fd_in);
  }
  if (fd_out != -1) {
    close(fd_out);
  }
  if (fd_err != -1) {
    close(fd_err);
  }
  int status = 0;
  if (waitpid(pid, &status, 0) == -1) {
    std::cerr << "RunExternalProgram failed in waitpid\n";
    std::cerr << "program=" << program << "\n";
    throw TerminalException{1};
  }
  if (WIFEXITED(status)) {
    return WEXITSTATUS(status);
  }
  if (WIFSIGNALED(status)) {
    return 128 + WTERMSIG(status);
  }
  return status;
}

// clang-format off
#endif  // SRC_BASIC_BASIC_EXTERNAL_PROGRAM_POSIX_H_
// clang-format on
