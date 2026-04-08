// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_FILE_H_
#define SRC_BASIC_BASIC_FILE_H_

#include "Basic_string.h"
#include "Temp_common.h"
#include <dirent.h>
#include <errno.h>
#include <filesystem>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
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

int RunExternalProgram(std::string const &program,
                       std::vector<std::string> const &inputs,
                       std::optional<std::string> const &output,
                       std::optional<std::string> const &error) {
  if (program.empty()) {
    std::cerr << "RunExternalProgram requires a non-empty program name\n";
    throw TerminalException{1};
  }
  int fd_out = -1;
  int fd_err = -1;
  if (output) {
    fd_out = open(output->c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_out == -1) {
      std::cerr << "RunExternalProgram failed to open output file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "output=" << *output << "\n";
      throw TerminalException{1};
    }
  }
  if (error) {
    fd_err = open(error->c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd_err == -1) {
      if (fd_out != -1) {
        close(fd_out);
      }
      std::cerr << "RunExternalProgram failed to open error file\n";
      std::cerr << "program=" << program << "\n";
      std::cerr << "error=" << *error << "\n";
      throw TerminalException{1};
    }
  }
  pid_t pid = fork();
  if (pid == -1) {
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
  std::string ePath = eDir + ".";
  DIR *dirp = opendir(ePath.c_str());
  if (dirp == nullptr) {
    std::cerr << "Error in routine FILE_GetDirectoryListFile\n";
    std::cerr << "Error in call to opendir\n";
    std::cerr << "eDir = " << eDir << "\n";
    throw TerminalException{1};
  }
  struct dirent *dp;
  std::vector<std::string> ListFile;
  while ((dp = readdir(dirp)) != nullptr) {
    std::string eName = dp->d_name;
    //    free(dp); // not sure this is portable
    if (eName != ".." && eName != ".")
      ListFile.push_back(eName);
  }
  int err = closedir(dirp);
  if (err != 0) {
    std::cerr << "err=" << err << "\n";
    printf("Oh dear, something went wrong with ls! %s\n", strerror(errno));
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
  int status;
  struct stat st_buf;
  status = stat(eFile.c_str(), &st_buf);
  if (status != 0) {
    std::cerr << "Problem in FILE_IsRegularFile\n";
    std::cerr << "Error, errno = " << errno << "\n";
    std::cerr << "eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  if (S_ISREG(st_buf.st_mode)) {
    return true;
  }
  return false;
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
  if (ProgName.find('/') != std::string::npos) {
    return access(ProgName.c_str(), X_OK) == 0;
  }
  char const *path_env = std::getenv("PATH");
  if (path_env == nullptr) {
    return false;
  }
  std::string path_str(path_env);
  size_t pos = 0;
  while (true) {
    size_t next = path_str.find(':', pos);
    std::string dir;
    if (next == std::string::npos) {
      dir = path_str.substr(pos);
    } else {
      dir = path_str.substr(pos, next - pos);
    }
    if (dir.empty()) {
      dir = ".";
    }
    std::filesystem::path candidate = std::filesystem::path(dir) / ProgName;
    if (access(candidate.c_str(), X_OK) == 0) {
      return true;
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
