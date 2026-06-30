// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_FILE_H_
#define SRC_BASIC_BASIC_FILE_H_

#include "Basic_string.h"
#include "Temp_common.h"
#include <filesystem>
#include <string>
#include <unistd.h>
#include <vector>

template <typename F>
void print_stderr_stdout_file(std::string const &FileOut, F f) {
  if (FileOut == "stderr")
    return f(std::cerr);
  if (FileOut == "stdout")
    return f(std::cout);
  std::ofstream os(FileOut);
  return f(os);
}

void CopyOperation(std::string const &SrcFile, std::string const &DstFile) {
  std::error_code ec;
  std::filesystem::copy_file(SrcFile, DstFile,
                             std::filesystem::copy_options::overwrite_existing,
                             ec);
  if (ec) {
    std::cerr << "Error in copy operation\n";
    std::cerr << "SrcFile=" << SrcFile << "\n";
    std::cerr << "DstFile=" << DstFile << "\n";
    std::cerr << "ec.message()=" << ec.message() << "\n";
    throw TerminalException{1};
  }
}

bool IsExistingFile(std::string const &eFile) {
  return std::ifstream(eFile).good();
}

std::string FindAvailableFileFromPrefix(std::string const &prefix) {
  size_t iFile = 0;
  while (true) {
    std::string FullFile = prefix + std::to_string(iFile);
    if (!IsExistingFile(FullFile)) {
      return FullFile;
    }
    iFile++;
  }
}

void IsExistingFileDie(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "Is missing. DIE\n";
    throw TerminalException{1};
  }
}

std::vector<std::string> ReadFullFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "ReadFullFile eFile=" << eFile << "\n";
    std::cerr << "Missing file\n";
    throw TerminalException{1};
  }
  std::ifstream is(eFile);
  std::string line;
  std::vector<std::string> ListLines;
  while (std::getline(is, line))
    ListLines.push_back(line);
  return ListLines;
}

std::string FILE_GetNakedFilename(std::string const &eFileFull) {
  size_t pos = eFileFull.rfind('/');
  if (pos == std::string::npos)
    return eFileFull;
  return eFileFull.substr(pos + 1);
}

std::string FILE_RemoveExtension(std::string const &eFileIn) {
  size_t pos = eFileIn.rfind('.');
  if (pos == std::string::npos) {
    std::cerr << "Failed to find the dot .\n";
    std::cerr << "eFileIn=" << eFileIn << "\n";
    throw TerminalException{1};
  }
  return eFileIn.substr(0, pos);
}

std::string FILE_GetDirectoryOfFileName(std::string const &eFileFull) {
  size_t pos = eFileFull.rfind('/');
  if (pos == std::string::npos) {
    std::cerr << "The file has no / so cannot find the directory\n";
    throw TerminalException{1};
  }
  return eFileFull.substr(0, pos + 1);
}

std::vector<std::string> FILE_GetDirectoryListFile(std::string const &eDir) {
  std::error_code ec;
  std::vector<std::string> ListFile;
  for (auto const &entry : std::filesystem::directory_iterator(eDir, ec)) {
    ListFile.push_back(entry.path().filename().string());
  }
  if (ec) {
    std::cerr << "Error in FILE_GetDirectoryListFile: " << ec.message() << "\n";
    std::cerr << "eDir = " << eDir << "\n";
    throw TerminalException{1};
  }
  return ListFile;
}

bool FILE_IsDirectoryEmpty(std::string const &eDir) {
  return FILE_GetDirectoryListFile(eDir).empty();
}

bool FILE_CheckFinalShashDirectory(std::string const &eDir) {
  return eDir.ends_with('/');
}

bool FILE_IsRegularFile(std::string const &eFile) {
  std::error_code ec;
  bool reg = std::filesystem::is_regular_file(eFile, ec);
  if (ec) {
    std::cerr << "Problem in FILE_IsRegularFile: " << ec.message() << "\n";
    std::cerr << "eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  return reg;
}

std::vector<std::string>
FILE_GetDirectoryFilesRecursively(std::string const &eDir) {
  std::vector<std::string> ListDir{eDir};
  std::vector<std::string> ListFile;
  while (true) {
    std::vector<std::string> NewListDir;
    for (auto &fDir : ListDir) {
      std::vector<std::string> LocalListFile = FILE_GetDirectoryListFile(fDir);
      for (auto &eFile : LocalListFile) {
        std::string NewEnt = fDir + eFile;
        if (FILE_IsRegularFile(NewEnt)) {
          ListFile.push_back(NewEnt);
        } else {
          std::string NewDir = NewEnt + "/";
          NewListDir.push_back(NewDir);
        }
      }
    }
    if (NewListDir.empty())
      break;
    ListDir = NewListDir;
  }
  return ListFile;
}

std::vector<std::string>
FILE_DirectoryFilesSpecificExtension(std::string const &ePrefix,
                                     std::string const &eExtension) {
  std::vector<std::string> ListFile =
      FILE_GetDirectoryFilesRecursively(ePrefix);
  std::vector<std::string> RetListFile;
  for (auto &eFile : ListFile) {
    if (eFile.ends_with(eExtension))
      RetListFile.push_back(eFile);
  }
  return RetListFile;
}

std::vector<std::string>
FILE_DirectoryMatchingPrefixExtension(std::string const &ePrefix,
                                      std::string const &eExtension) {
  std::string eDir = FILE_GetDirectoryOfFileName(ePrefix);
  std::string eBeginStr = FILE_GetNakedFilename(ePrefix);
  std::cerr << "ePrefix=" << ePrefix << "\n";
  std::cerr << "eDir=" << eDir << " eBeginStr=" << eBeginStr << "\n";
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(eDir);
  std::vector<std::string> ListFile_RET;
  for (auto &eFile : ListFile) {
    if (eFile.starts_with(eBeginStr) && eFile.ends_with(eExtension)) {
      ListFile_RET.push_back(eDir + eFile);
    }
  }
  return ListFile_RET;
}

// If the ePrefix ends with / then we do recursive listing
// If the output is of the form path/WWM_output_
// then returns all the files having the
std::vector<std::string>
FILE_DirectoryFilesSpecificExtension_Gen(std::string const &ePrefix,
                                         std::string const &eExtension) {
  size_t len = ePrefix.size();
  if (len == 0) {
    std::cerr << "FILE_DirectoryFilesSpecificExtension_Gen requires a non-empty "
                 "prefix\n";
    throw TerminalException{1};
  }
  char LastChar = ePrefix[len - 1];
  if (LastChar == '/')
    return FILE_DirectoryFilesSpecificExtension(ePrefix, eExtension);
  //
  return FILE_DirectoryMatchingPrefixExtension(ePrefix, eExtension);
}

bool IsExistingDirectory(std::string const &ThePrefix) {
  if (0 != access(ThePrefix.c_str(), F_OK)) {
    if (ENOENT == errno) {
      // does not exist
      return false;
    }
    if (ENOTDIR == errno) {
      return false;
      // not a directory
    }
    std::cerr << "Should not happen a priori\n";
    throw TerminalException{1};
  }
  return true;
}

void RemoveEmptyDirectory(std::string const &eDir) {
  std::error_code ec;
  bool removed = std::filesystem::remove(eDir, ec);
  if (ec || !removed) {
    std::cerr << "Error in RemoveEmptyDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    if (ec) {
      std::cerr << "ec.message()=" << ec.message() << "\n";
    } else {
      std::cerr << "Directory was not removed\n";
    }
    throw TerminalException{1};
  }
}

void RemoveFile(std::string const &eFile) { std::remove(eFile.c_str()); }

void RemoveFileIfExist(std::string const &eFile) {
  if (IsExistingFile(eFile))
    RemoveFile(eFile);
}

bool IsProgramInPath(std::string const &ProgName) {
  if (ProgName.empty()) {
    return false;
  }
  namespace fs = std::filesystem;
  // Platform facts derived from std::filesystem rather than preprocessor:
  // path::preferred_separator is '/' on POSIX-style systems and L'\\' on
  // Windows-style ones. From that we pick the PATH list separator and the
  // set of executable extensions to probe.
  constexpr bool windows_style = fs::path::preferred_separator != '/';
  constexpr char path_list_separator = windows_style ? ';' : ':';
  static const std::vector<std::string> exec_exts =
      windows_style
          ? std::vector<std::string>{"", ".exe", ".com", ".bat", ".cmd"}
          : std::vector<std::string>{""};
  auto is_runnable = [](fs::path const &p) -> bool {
    std::error_code ec;
    fs::file_status st = fs::status(p, ec);
    if (ec || !fs::is_regular_file(st)) {
      return false;
    }
    // On POSIX this checks the executable bit; on Windows libstdc++ reports
    // owner_exec for every regular file, which collapses to an existence
    // check (executability is implied by the extension list above).
    return (st.permissions() & fs::perms::owner_exec) != fs::perms::none;
  };
  if (fs::path(ProgName).has_parent_path()) {
    for (auto const &ext : exec_exts) {
      if (is_runnable(fs::path(ProgName + ext))) {
        return true;
      }
    }
    return false;
  }
  char const *path_env = std::getenv("PATH");
  if (path_env == nullptr) {
    return false;
  }
  std::string path_str(path_env);
  size_t pos = 0;
  while (true) {
    size_t next = path_str.find(path_list_separator, pos);
    std::string dir;
    if (next == std::string::npos) {
      dir = path_str.substr(pos);
    } else {
      dir = path_str.substr(pos, next - pos);
    }
    if (dir.empty()) {
      dir = ".";
    }
    for (auto const &ext : exec_exts) {
      fs::path candidate = fs::path(dir) / (ProgName + ext);
      if (is_runnable(candidate)) {
        return true;
      }
    }
    if (next == std::string::npos) {
      break;
    }
    pos = next + 1;
  }
  return false;
}

std::string FILE_RemoveEndingExtension(std::string const &FileName,
                                       std::string const &TheExtension) {
  size_t len = FileName.size();
  std::ptrdiff_t iCharLast = -1;
  for (size_t iChar = 0; iChar < len; iChar++) {
    char eChar = FileName[iChar];
    if (eChar == '.')
      iCharLast = std::ptrdiff_t(iChar);
  }
  if (iCharLast == -1)
    return FileName;
  std::string FileNameRed = FileName.substr(0, size_t(iCharLast));
  std::string eExtension =
      FileName.substr(size_t(iCharLast + 1), size_t(len - 1 - iCharLast));
  if (eExtension == TheExtension)
    return FileNameRed;
  return FileName;
}

void RemoveFileSpecificExtension(std::string const &ThePrefix,
                                 std::string const &TheExtension) {
  bool test = IsExistingDirectory(ThePrefix);
  if (!test)
    return;
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(ThePrefix);
  for (auto &eFile : ListFile) {
    if (eFile.ends_with(TheExtension)) {
      std::string eFileTot = ThePrefix + eFile;
      RemoveFile(eFileTot);
    }
  }
}

void RemoveFileInDirectory(std::string const &ThePrefix) {
  bool test = IsExistingDirectory(ThePrefix);
  if (!test)
    return;
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(ThePrefix);
  for (auto &eFile : ListFile) {
    std::string eFileTot = ThePrefix + eFile;
    if (!FILE_IsRegularFile(eFileTot)) {
      std::cerr << "We have subdirectories in ThePrefix=" << ThePrefix << "\n";
      std::cerr
          << "Therefore it seems unwise to remove automatically by program\n";
      std::cerr << "the entries in the directory since it is far too likely to "
                   "destroy valuable file.\n";
      std::cerr << "If not an error, please remove things manually\n";
      throw TerminalException{1};
    }
  }
  for (auto &eFile : ListFile) {
    std::string eFileTot = ThePrefix + eFile;
    RemoveFile(eFileTot);
  }
}

int FILE_GetNumberLine(std::string const &eFile) {
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(eFile);
  while (std::getline(myfile, line))
    ++number_of_lines;
  return number_of_lines;
}

std::string GetCurrentDirectory() {
  std::error_code ec;
  std::filesystem::path ePath = std::filesystem::current_path(ec);
  if (ec) {
    std::cerr << "Error while trying to use current_path\n";
    std::cerr << "ec.message()=" << ec.message() << "\n";
    throw TerminalException{1};
  }
  std::string eRet = ePath.string();
  eRet += "/";
  return eRet;
}

std::string FILE_GetAbsoluteDirectory(std::string const &ePrefix) {
  char FirstChar = ePrefix[0];
  if (FirstChar == '/') {
    return ePrefix;
  } else {
    std::string ePWD = GetCurrentDirectory();
    return ePWD + ePrefix;
  }
}

bool FILE_CheckPrefix(std::string const &ePrefix) {
  return ePrefix.ends_with('/');
}

std::string ExtractDirectoryFromFileString(std::string const &eFile) {
  size_t pos = eFile.rfind('/');
  if (pos == std::string::npos) {
    std::cerr << "Error in ExtractDirectoryFromFileString\n";
    throw TerminalException{1};
  }
  return eFile.substr(0, pos + 1);
}

bool FILE_IsFileMakeable(std::string const &eFile) {
  std::string eDir = ExtractDirectoryFromFileString(eFile);
  if (!IsExistingFile(eDir))
    return false;
  return true;
}

void CreateDirectory(std::string const &eDir) {
  if (eDir.empty()) {
    std::cerr << "Error in CreateDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    throw TerminalException{1};
  }
  std::error_code ec;
  std::filesystem::create_directories(eDir, ec);
  if (ec && !std::filesystem::is_directory(eDir)) {
    std::cerr << "Error in CreateDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    std::cerr << "ec.message()=" << ec.message() << "\n";
    throw TerminalException{1};
  }
}

std::vector<std::string> ls_operation(std::string const &ThePrefix) {
  std::error_code ec;
  std::vector<std::string> ListFile;
  for (auto const &entry : std::filesystem::directory_iterator(ThePrefix, ec)) {
    if (ec) {
      break;
    }
    std::string eFile = entry.path().filename().string();
    if (!eFile.empty()) {
      ListFile.push_back(eFile);
    }
  }
  if (ec) {
    std::cerr << "Error in ls_operation\n";
    std::cerr << "ThePrefix=" << ThePrefix << "\n";
    std::cerr << "ec.message()=" << ec.message() << "\n";
    throw TerminalException{1};
  }
  return ListFile;
}

struct TempFile {
private:
  std::string FileName;

public:
  TempFile() = delete;
  TempFile(std::string const &eFile) { FileName = eFile; }
  TempFile(char *eFile) { FileName = eFile; }
  ~TempFile() {
    if (IsExistingFile(FileName)) {
      RemoveFile(FileName);
    }
  }
  //
  bool IsExisting() const { return IsExistingFile(FileName); }
  std::string string() const { return FileName; }
};

// This does not provide all the facility that we are after
// The problem is that when we have an exit(1)
// the destructors of the variables are not called.
//
// This means that we need to put a throw statement.
// in order to have the temporary directory eliminated.
// So eliminate all the exit(1) statement and replace by throw
// and catch the throws.
struct CondTempDirectory {
private:
  bool used;
  std::string DirName;

public:
  CondTempDirectory() {
    used = false;
    DirName = "unset_and_irrelevant";
  }
  CondTempDirectory(bool const &eUsed, std::string const &eDir) {
    used = eUsed;
    if (used) {
      DirName = eDir;
      CreateDirectory(DirName);
    } else {
      DirName = "unset_and_irrelevant";
    }
  }
  CondTempDirectory &operator=(CondTempDirectory &&eTemp) {
    used = eTemp.usedness();
    DirName = eTemp.str();
    eTemp.DirName = "unset_and_irrelevant";
    return *this;
  }
  CondTempDirectory(CondTempDirectory &&eTemp)
      : used(eTemp.usedness()), DirName(eTemp.str()) {}
  ~CondTempDirectory() {
    if (used && DirName != "unset_and_irrelevant") {
      if (IsExistingDirectory(DirName)) {
        if (!FILE_IsDirectoryEmpty(DirName)) {
          std::cerr << "Keeping " << DirName << " since it is not empty\n";
        } else {
          RemoveFile(DirName);
        }
      }
    }
  }
  //
  bool IsExisting() const { return IsExistingDirectory(DirName); }
  bool usedness() const { return used; }
  std::string str() const { return DirName; }
};

// clang-format off
#endif  // SRC_BASIC_BASIC_FILE_H_
// clang-format on
