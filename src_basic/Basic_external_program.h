// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_EXTERNAL_PROGRAM_H_
#define SRC_BASIC_BASIC_EXTERNAL_PROGRAM_H_

// Runs an external program in a child process and returns its exit code:
//
//   int RunExternalProgram(std::string const &program,
//                          std::vector<std::string> const &inputs,
//                          std::optional<std::string> const &std_input,
//                          std::optional<std::string> const &std_output,
//                          std::optional<std::string> const &std_error);
//
// The implementation is platform specific (POSIX fork / execvp / waitpid versus
// the Windows CRT _spawnvp facility), so select the right one here. Include this
// header, never the _posix / _windows variants directly.

#if defined(_WIN32)
#include "Basic_external_program_windows.h"
#else
#include "Basic_external_program_posix.h"
#endif

// clang-format off
#endif  // SRC_BASIC_BASIC_EXTERNAL_PROGRAM_H_
// clang-format on
