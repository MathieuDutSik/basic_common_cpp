// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Basic_file.h"
#include <algorithm>
#include <chrono>
#include <filesystem>

void check(bool cond, std::string const &msg) {
  if (!cond) {
    std::cerr << "Filesystem test failure: " << msg << "\n";
    throw TerminalException{1};
  }
}

int main() {
  try {
    namespace fs = std::filesystem;

    std::string current_dir = GetCurrentDirectory();
    std::string expected_current_dir = fs::current_path().string() + "/";
    check(current_dir == expected_current_dir,
          "GetCurrentDirectory should match std::filesystem::current_path");

    auto now = std::chrono::steady_clock::now().time_since_epoch().count();
    fs::path root =
        fs::temp_directory_path() / ("basic_common_cpp_fs_test_" + std::to_string(now));
    fs::remove_all(root);

    std::string root_dir = root.string() + "/";
    std::string nested_dir = root_dir + "nested/inner/";
    std::string flat_dir = root_dir + "flat/";

    CreateDirectory(nested_dir);
    CreateDirectory(flat_dir);

    check(IsExistingDirectory(root_dir), "root directory should exist");
    check(IsExistingDirectory(nested_dir), "nested directory should exist");
    check(FILE_IsDirectoryEmpty(nested_dir), "nested directory should start empty");

    std::string source_file = nested_dir + "source.txt";
    {
      std::ofstream os(source_file);
      os << "alpha\n";
      os << "beta\n";
    }

    check(IsExistingFile(source_file), "source file should exist");
    check(FILE_IsRegularFile(source_file), "source should be a regular file");
    check(!FILE_IsDirectoryEmpty(nested_dir), "nested directory should no longer be empty");

    std::vector<std::string> lines = ReadFullFile(source_file);
    check(lines.size() == 2, "source file should have two lines");
    check(lines[0] == "alpha" && lines[1] == "beta",
          "source file contents should match");

    std::string copied_file = flat_dir + "copied.txt";
    CopyOperation(source_file, copied_file);
    check(IsExistingFile(copied_file), "copied file should exist");
    check(ReadFullFile(copied_file) == lines, "copied file contents should match");

    std::string matching_file = flat_dir + "sample.log";
    {
      std::ofstream os(matching_file);
      os << "log\n";
    }

    std::vector<std::string> top_entries = FILE_GetDirectoryListFile(root_dir);
    std::sort(top_entries.begin(), top_entries.end());
    check(top_entries.size() == 2, "root directory should contain two entries");
    check(top_entries[0] == "flat" && top_entries[1] == "nested",
          "root directory entries should match");

    std::vector<std::string> ls_entries = ls_operation(root_dir);
    std::sort(ls_entries.begin(), ls_entries.end());
    check(ls_entries.size() == 2, "ls_operation should find two root entries");
    check(ls_entries[0] == "flat" && ls_entries[1] == "nested",
          "ls_operation entries should match");

    std::vector<std::string> recursive_files =
        FILE_GetDirectoryFilesRecursively(root_dir);
    std::sort(recursive_files.begin(), recursive_files.end());
    check(recursive_files.size() == 3,
          "recursive listing should find three files");
    check(recursive_files[0] == flat_dir + "copied.txt" &&
              recursive_files[1] == flat_dir + "sample.log" &&
              recursive_files[2] == nested_dir + "source.txt",
          "recursive listing should return the expected files");

    std::vector<std::string> txt_files =
        FILE_DirectoryFilesSpecificExtension_Gen(root_dir, "txt");
    std::sort(txt_files.begin(), txt_files.end());
    check(txt_files.size() == 2, "txt listing should find two files");
    check(txt_files[0] == flat_dir + "copied.txt" &&
              txt_files[1] == nested_dir + "source.txt",
          "txt listing should match expected files");

    std::vector<std::string> prefix_files =
        FILE_DirectoryFilesSpecificExtension_Gen(flat_dir + "cop", "txt");
    check(prefix_files.size() == 1 && prefix_files[0] == copied_file,
          "prefix listing should find copied.txt");

    RemoveFileSpecificExtension(flat_dir, "log");
    check(!IsExistingFile(matching_file), "log file should be removed");
    check(IsExistingFile(copied_file), "txt file should remain after log removal");

    RemoveFileInDirectory(flat_dir);
    check(FILE_IsDirectoryEmpty(flat_dir), "flat directory should be empty after removal");

    RemoveFile(source_file);
    check(FILE_IsDirectoryEmpty(nested_dir), "nested directory should be empty after file removal");

    RemoveEmptyDirectory(nested_dir);
    check(!IsExistingDirectory(nested_dir), "nested directory should be removed");

    RemoveEmptyDirectory(root_dir + "nested/");
    RemoveEmptyDirectory(flat_dir);
    RemoveEmptyDirectory(root_dir);

    check(!IsExistingDirectory(root_dir), "root directory should be removed");
    std::cerr << "Normal termination of Test_filesystem\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_filesystem\n";
    exit(e.eVal);
  }
}
