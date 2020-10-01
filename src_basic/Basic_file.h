#ifndef BASIC_FILE_INCLUDE
#define BASIC_FILE_INCLUDE

#include "Temp_common.h"
#include "Basic_string.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


void CopyOperation(std::string const& SrcFile, std::string const& DstFile)
{
  std::string eComm="cp " + SrcFile + " " + DstFile;
  int iret=system(eComm.c_str() );
  if (iret != 0) {
    std::cerr << "Error in copy operation\n";
    std::cerr << "SrcFile=" << SrcFile << "\n";
    std::cerr << "DstFile=" << DstFile << "\n";
    throw TerminalException{1};
  }
}


bool IsExistingFile(std::string const& eFile)
{
  std::ifstream f(eFile.c_str());
  if (f.good()) {
    f.close();
    return true;
  } else {
    f.close();
    return false;
  }   
}

void IsExistingFileDie(std::string const& eFile)
{
  if (!IsExistingFile(eFile)) {
    std::cerr << "The file eFile = " << eFile << "\n";
    std::cerr << "Is missing. DIE\n";
    throw TerminalException{1};
  }
}





std::vector<std::string> ReadFullFile(std::string const& eFile)
{
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



std::string FILE_GetNakedFilename(std::string const& eFileFull)
{
  std::vector<std::string> LStr=STRING_Split(eFileFull, "/");
  std::string eLast=LStr[LStr.size() -1];
  return eLast;
}

std::string FILE_RemoveExtension(std::string const& eFileIn)
{
  std::ptrdiff_t pos=-1;
  size_t len=eFileIn.size();
  std::string eChar=".";
  for (size_t i=0; i<len; i++) {
    if (eFileIn.substr(i,1) == eChar) {
      pos=i;
    }
  }
  if (pos == -1) {
    std::cerr << "Failed to find the dot .\n";
    std::cerr << "eFileIn=" << eFileIn << "\n";
    throw TerminalException{1};
  }
  return eFileIn.substr(0,size_t(pos));
}



std::string FILE_GetDirectoryOfFileName(std::string const& eFileFull)
{
  size_t len=eFileFull.size();
  std::ptrdiff_t posfound=-1;
  for (size_t u=0; u<len; u++) {
    std::string eChar=eFileFull.substr(u,1);
    if (eChar == "/")
      posfound=std::ptrdiff_t(u);
  }
  if (posfound == -1) {
    std::cerr << "The file has no / so cannot find the directory\n";
    throw TerminalException{1};
  }
  return eFileFull.substr(0, size_t(posfound+1));
}





std::vector<std::string> FILE_GetDirectoryListFile(std::string const& eDir)
{
  std::string ePath=eDir + ".";
  DIR* dirp=opendir(ePath.c_str());
  if (dirp == NULL) {
    std::cerr << "Error in routine FILE_GetDirectoryListFile\n";
    std::cerr << "Error in call to opendir\n";
    std::cerr << "eDir = " << eDir << "\n";
    throw TerminalException{1};
  }
  struct dirent *dp;
  std::vector<std::string> ListFile;
  while ((dp = readdir(dirp)) != NULL) {
    std::string eName=dp->d_name;
    //    free(dp); // not sure this is portable
    if (eName != ".." && eName != ".")
      ListFile.push_back(eName);
  }
  int err=closedir(dirp);
  if (err != 0) {
    std::cerr << "err=" << err << "\n";
    printf("Oh dear, something went wrong with ls! %s\n", strerror(errno));
    throw TerminalException{1};
  }
  return ListFile;
}

bool FILE_IsDirectoryEmpty(std::string const& eDir)
{
  std::vector<std::string> TheList = FILE_GetDirectoryListFile(eDir);
  if (TheList.size() == 0)
    return true;
  return false;
}



bool FILE_CheckFinalShashDirectory(std::string const& eDir)
{
  size_t len=eDir.size();
  std::string eChar=eDir.substr(len-1,1);
  if (eChar == "/")
    return true;
  return false;
}



bool FILE_IsRegularFile(std::string const& eFile)
{
  int status;
  struct stat st_buf;  
  status = stat(eFile.c_str(), &st_buf);
  if (status != 0) {
    std::cerr << "Problem in FILE_IsRegularFile\n";
    std::cerr << "Error, errno = " << errno << "\n";
    std::cerr << "eFile=" << eFile << "\n";
    throw TerminalException{1};
  }
  if (S_ISREG (st_buf.st_mode)) {
    return true;
  }
  return false;
}



std::vector<std::string> FILE_GetDirectoryFilesRecursively(std::string const& eDir)
{
  //  std::cerr << "Beginning of FILE_GetDirectoryFilesRecursively\n";
  std::vector<std::string> ListDir{eDir};
  std::vector<std::string> ListFile;
  while(true) {
    std::vector<std::string> NewListDir;
    for (auto & fDir : ListDir) {
      std::vector<std::string> LocalListFile=FILE_GetDirectoryListFile(fDir);
      for (auto & eFile : LocalListFile) {
	std::string NewEnt=fDir + eFile;
	if (FILE_IsRegularFile(NewEnt)) {
	  ListFile.push_back(NewEnt);
	}
	else {
	  std::string NewDir=NewEnt + "/";
	  NewListDir.push_back(NewDir);
	}
      }
    }
    if (NewListDir.size() == 0)
      break;
    ListDir=NewListDir;
  }
  //  std::cerr << "Ending of FILE_GetDirectoryFilesRecursively\n";
  return ListFile;
}

std::vector<std::string> FILE_DirectoryFilesSpecificExtension(std::string const& ePrefix, std::string const& eExtension)
{
  std::vector<std::string> ListFile=FILE_GetDirectoryFilesRecursively(ePrefix);
  std::vector<std::string> RetListFile;
  size_t lenExt=eExtension.size();
  for (auto & eFile : ListFile) {
    size_t len=eFile.size();
    if (len > lenExt) {
      std::string eEnd=eFile.substr(len-lenExt,lenExt);
      if (eEnd == eExtension)
	RetListFile.push_back(eFile);
    }
  }
  return RetListFile;
}



std::vector<std::string> FILE_DirectoryMatchingPrefixExtension(std::string const& ePrefix, std::string const& eExtension)
{
  std::string eDir = FILE_GetDirectoryOfFileName(ePrefix);
  std::string eBeginStr = FILE_GetNakedFilename(ePrefix);
  std::cerr << "ePrefix=" << ePrefix << "\n";
  std::cerr << "eDir=" << eDir << " eBeginStr=" << eBeginStr << "\n";
  size_t lenExt=eExtension.size();
  size_t lenBegin=eBeginStr.size();
  std::vector<std::string> ListFile = FILE_GetDirectoryListFile(eDir);
  std::vector<std::string> ListFile_RET;
  for (auto & eFile : ListFile) {
    size_t len = eFile.size();
    if (len > lenBegin && len > lenExt) {
      std::string str1=eFile.substr(0, lenBegin);
      std::string str2=eFile.substr(len-lenExt, lenExt);
      if (str1 == eBeginStr && str2 == eExtension) {
	std::string NewFile = eDir + eFile;
	ListFile_RET.push_back(NewFile);
      }
    }
  }
  return ListFile_RET;
}



// If the ePrefix ends with / then we do recursive listing
// If the output is of the form path/WWM_output_
// then returns all the files having the 
std::vector<std::string> FILE_DirectoryFilesSpecificExtension_Gen(std::string const& ePrefix, std::string const& eExtension)
{
  size_t len=ePrefix.size();
  std::string LastChar = ePrefix.substr(len-1,1);
  if (LastChar == std::string("/"))
    return FILE_DirectoryFilesSpecificExtension(ePrefix, eExtension);
  //
  return FILE_DirectoryMatchingPrefixExtension(ePrefix, eExtension);
}





bool IsExistingDirectory(std::string const& ThePrefix)
{
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


void RemoveEmptyDirectory(std::string const& eDir)
{
  std::string eOrder="rm -d " + eDir;
  int iret=system(eOrder.c_str() );
  if (iret == -1) {
    std::cerr << "Error in RemoveEmptyDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    std::cerr << "eOrder=" << eOrder << "\n";
    std::cerr << "unable to run the process\n";
    throw TerminalException{1};
  }
}




void RemoveFile(std::string const& eFile)
{
  std::remove(eFile.c_str());
}

void RemoveFileIfExist(std::string const& eFile)
{
  if (IsExistingFile(eFile))
    RemoveFile(eFile);
}


std::string FILE_RemoveEndingExtension(std::string const& FileName, std::string const& TheExtension)
{
  size_t len=FileName.size();
  std::ptrdiff_t iCharLast=-1;
  for (size_t iChar=0; iChar<len; iChar++) {
    std::string eChar=FileName.substr(iChar,1);
    if (eChar == ".")
      iCharLast=std::ptrdiff_t(iChar);
  }
  if (iCharLast == -1)
    return FileName;
  std::string FileNameRed=FileName.substr(0,size_t(iCharLast));
  std::string eExtension=FileName.substr(size_t(iCharLast+1),size_t(len-1-iCharLast));
  if (eExtension == TheExtension)
    return FileNameRed;
  return FileName;
}


void RemoveFileSpecificExtension(std::string const& ThePrefix, std::string const& TheExtension)
{
  bool test=IsExistingDirectory(ThePrefix);
  if (!test)
    return;
  std::vector<std::string> ListFile=FILE_GetDirectoryListFile(ThePrefix);
  size_t nbCharEnd=TheExtension.size();
  for (auto & eFile : ListFile) {
    size_t len=eFile.size();
    if (len > nbCharEnd) {
      std::string TheEnd=eFile.substr(len-nbCharEnd, nbCharEnd);
      if (TheEnd == TheExtension) {
	std::string eFileTot=ThePrefix + eFile;
	RemoveFile(eFileTot);
      }
    }
  }
}


void RemoveFileInDirectory(std::string const& ThePrefix)
{
  bool test=IsExistingDirectory(ThePrefix);
  if (!test)
    return;
  std::vector<std::string> ListFile=FILE_GetDirectoryListFile(ThePrefix);
  for (auto & eFile : ListFile) {
    std::string eFileTot=ThePrefix + eFile;
    if (!FILE_IsRegularFile(eFileTot)) {
      std::cerr << "We have subdirectories in ThePrefix=" << ThePrefix << "\n";
      std::cerr << "Therefore it seems unwise to remove automatically by program\n";
      std::cerr << "the entries in the directory since it is far too likely to destroy valuable file.\n";
      std::cerr << "If not an error, please remove things manually\n";
      throw TerminalException{1};
    }
  }
  for (auto & eFile : ListFile) {
    std::string eFileTot=ThePrefix + eFile;
    RemoveFile(eFileTot);
  }
}





int FILE_GetNumberLine(std::string const& eFile)
{
  int number_of_lines = 0;
  std::string line;
  std::ifstream myfile(eFile);
  while (std::getline(myfile, line))
    ++number_of_lines;
  return number_of_lines;
}




#ifndef WINDOWS
std::string GetCurrentDirectory()
{
  int size = pathconf(".", _PC_PATH_MAX);
  std::vector<char> buf(size);
  char* ptr=getcwd(buf.data(), (size_t)size);
  if (ptr == NULL && errno != ERANGE) {
    std::cerr << "Error while trying to use getcwd\n";
    throw TerminalException{1};
  }
  std::string eRet = buf.data();
  eRet=eRet + "/";
  if (ptr != NULL) {
    if (ptr != buf.data()) {
      std::cerr << "Before ptr freeing\n";
      std::free(ptr);
      std::cerr << "After ptr freeing\n";
    }
  }
  return eRet;
}
#endif


#ifndef WINDOWS
std::string FILE_GetAbsoluteDirectory(std::string const& ePrefix)
{
  std::string FirstChar=ePrefix.substr(0, 1);
  if (FirstChar == "/") {
    return ePrefix;
  }
  else {
    std::string ePWD=GetCurrentDirectory();
    return ePWD + ePrefix;
  }
}
#endif

bool FILE_CheckPrefix(std::string const& ePrefix)
{
  size_t len=ePrefix.size();
  std::string LastChar=ePrefix.substr(len-1,1);
  //  std::cerr << "LastChar=" << LastChar << "\n";
  if (LastChar == "/")
    return true;
  return false;
}

std::string ExtractDirectoryFromFileString(std::string const& eFile)
{
  size_t len=eFile.size();
  std::ptrdiff_t iCharFinal=-1;
  for (size_t iChar=0; iChar<len; iChar++) {
    std::string eChar=eFile.substr(iChar,1);
    if (eChar == "/")
      iCharFinal=std::ptrdiff_t(iChar);
  }
  if (iCharFinal == -1) {
    std::cerr << "Error in ExtractDirectoryFromFileString\n";
    throw TerminalException{1};
  }
  return eFile.substr(0,size_t(iCharFinal+1));
}


bool FILE_IsFileMakeable(std::string const& eFile)
{
  std::string eDir=ExtractDirectoryFromFileString(eFile);
  if (!IsExistingFile(eDir))
    return false;
  return true;
}





void CreateDirectory(std::string const& eDir)
{
  //  std::cerr << "eDir=" << eDir << "\n";
  const char *dir=eDir.c_str();
  char tmp[256];
  char *p = NULL;
  size_t len;
  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/')
    tmp[len - 1] = 0;
  for(p = tmp + 1; *p; p++)
    if(*p == '/') {
      *p = 0;
      //      std::cerr << "Before mkdir, tmp=" << tmp << "\n";
      mkdir(tmp, S_IRWXU);
      *p = '/';
    }
  //  std::cerr << "Before mkdir, tmp=" << tmp << "\n";
  mkdir(tmp, S_IRWXU);
}


void CreateDirectory_V1(std::string const& eDir)
{
  std::string eOrder="/bin/mkdir -p " + eDir;
  int iret=system(eOrder.c_str() );
  if (iret == -1) {
    std::cerr << "Error in CreateDirectory\n";
    std::cerr << "eDir=" << eDir << "\n";
    std::cerr << "eOrder=" << eOrder << "\n";
    std::cerr << "unable to run the process\n";
    throw TerminalException{1};
  }
}


std::vector<std::string> ls_operation(std::string const& ThePrefix)
{
  std::cerr << "Doing ls_operation\n";
  std::string strRand=random_string(20);
  std::string TmpFile="/tmp/file" + strRand;
  std::string ErrFile="/tmp/file" + strRand + ".err";
  std::string eOrder="ls " + ThePrefix + " > " + TmpFile + " 2> " + ErrFile;
  int iret=system(eOrder.c_str() );
  std::cerr << "iret=" << iret << "\n";
  if (iret != -1) {
    std::cerr << "Error in ls_operation\n";
    std::cerr << "ThePrefix=" << ThePrefix << "\n";
    std::cerr << "unable to run the process\n";
    throw TerminalException{1};
  }
  std::ifstream iserr(ErrFile);
  int nbCharErr=0;
  while(!iserr.eof()) {
    std::string PreStr;
    std::getline(iserr, PreStr);
    nbCharErr += PreStr.size();
  }
  if (nbCharErr > 0) {
    std::cerr << "Error in ls_operation\n";
    std::cerr << "We have nbCharErr = " << nbCharErr << "\n";
    std::cerr << "TmpFile=" << TmpFile << "\n";
    std::cerr << "ErrFile=" << ErrFile << "\n";
    std::cerr << "ThePrefix=" << ThePrefix << "\n";
    throw TerminalException{1};
  }
  //
  std::ifstream is(TmpFile);
  std::vector<std::string> ListFile;
  while(true) {
    if (is.eof()) {
      is.close();
      RemoveFile(TmpFile);
      RemoveFile(ErrFile);
      return ListFile;
    }
    std::string eFile;
    is >> eFile;
    ListFile.push_back(eFile);
  }
}



struct TempFile {
private:
  std::string FileName;
public:
  TempFile() = delete;
  TempFile(std::string const& eFile)
  {
    FileName=eFile;
  }
  TempFile(char* eFile)
  {
    FileName=eFile;
  }
  ~TempFile()
  {
    if (IsExistingFile(FileName)) {
      RemoveFile(FileName);
    }
  }
  //
  bool IsExisting() const
  {
    return IsExistingFile(FileName);
  }
  std::string string() const
  {
    return FileName;
  }
};



// This does not provide all the facility that we are after
// The problem is that when we have an exit(1)
// the destructors of the variables are not called.
//
// This means that we need to put a throw statement.
// in order to have the temporary directory eliminated.
// So eliminate all the exit(1) statement and replace by throw
// and catch the throws.
struct TempDirectory {
private:
  std::string DirName;
  bool IsInitialized;
public:
  TempDirectory()
  {
    DirName="unset_and_irrelevant";
    IsInitialized=false;
  }
  TempDirectory(std::string const& eDir)
  {
    DirName=eDir;
    CreateDirectory(DirName);
    IsInitialized=true;
  }
  TempDirectory& operator=(TempDirectory&& eTemp)
  {
    DirName=eTemp.str();
    IsInitialized=true;
    eTemp.DirName="unset_and_irrelevant";
    eTemp.IsInitialized=false;
    return *this;
  }
  TempDirectory(TempDirectory && eTemp) : DirName(eTemp.str()), IsInitialized(true)
  {
    eTemp.DirName="unset_and_irrelevant";
    eTemp.IsInitialized=false;
  }
  
  TempDirectory(const TempDirectory & eTemp) = delete;
  
  TempDirectory& operator=(const TempDirectory & eTemp) = delete;
  
  ~TempDirectory()
  {
    //    std::cerr << "Calling destructor\n";
    if (IsInitialized) {
      //      std::cerr << "  Destructor is really needed\n";
      //      std::cerr << "  DirName=" << DirName << "\n";
      if (IsExistingDirectory(DirName)) {
	//	std::cerr << "    Directory really exist\n";
	if (!FILE_IsDirectoryEmpty(DirName)) {
	  std::cerr << "Keeping " << DirName << " since it is not empty\n";
	}
	else {
	  RemoveFile(DirName);
	}
      }
    }
  }
  //
  bool IsExisting() const
  {
    return IsExistingDirectory(DirName);
  }
  std::string str() const
  {
    return DirName;
  }
};


struct CondTempDirectory {
private:
  bool used;
  std::string DirName;
public:
  CondTempDirectory()
  {
    used=false;
    DirName="unset_and_irrelevant";
  }
  CondTempDirectory(bool const& eUsed, std::string const& eDir)
  {
    used=eUsed;
    if (used) {
      DirName=eDir;
      CreateDirectory(DirName);
    }
    else {
      DirName="unset_and_irrelevant";
    }
  }
  CondTempDirectory& operator=(CondTempDirectory&& eTemp)
  {
    used=eTemp.usedness();
    DirName=eTemp.str();
    eTemp.DirName="unset_and_irrelevant";
    return *this;
  }
  CondTempDirectory(CondTempDirectory && eTemp) : used(eTemp.usedness()), DirName(eTemp.str())
  {
  }
  ~CondTempDirectory()
  {
    //    std::cerr << "Calling destructor\n";
    if (used && DirName != "unset_and_irrelevant") {
      //      std::cerr << "  Destructor is really needed\n";
      //      std::cerr << "  DirName=" << DirName << "\n";
      if (IsExistingDirectory(DirName)) {
	//	std::cerr << "    Directory really exist\n";
	if (!FILE_IsDirectoryEmpty(DirName)) {
	  std::cerr << "Keeping " << DirName << " since it is not empty\n";
	}
	else {
	  RemoveFile(DirName);
	}
      }
    }
  }
  //
  bool IsExisting() const
  {
    return IsExistingDirectory(DirName);
  }
  bool usedness() const
  {
    return used;
  }
  std::string str() const
  {
    return DirName;
  }
};



#endif
