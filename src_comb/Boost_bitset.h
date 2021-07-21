#ifndef INCLUDE_FACE_BITSET
#define INCLUDE_FACE_BITSET

// Boost libraries

#include "Temp_common.h"
#include "Boost_bitset_kernel.h"



std::vector<int> FaceToVector(Face const& eSet)
{
  size_t nbVert=eSet.count();
  std::vector<int> eList(nbVert);
  boost::dynamic_bitset<>::size_type aRow=eSet.find_first();
  for (size_t i=0; i<nbVert; i++) {
    eList[i]=int(aRow);
    aRow = eSet.find_next(aRow);
  }
  return eList;
}


std::vector<int> FaceTo01vector(Face const& eSet)
{
  size_t nbVert=eSet.size();
  size_t siz=eSet.count();
  std::vector<int> eList(nbVert,0);
  boost::dynamic_bitset<>::size_type aRow=eSet.find_first();
  for (size_t i=0; i<siz; i++) {
    eList[aRow]=1;
    aRow = eSet.find_next(aRow);
  }
  return eList;
}







void WriteFace(std::ostream & os, Face const& eList)
{
  size_t len = eList.size();
  os << len;
  for (size_t i=0; i<len; i++) {
    int eVal=eList[i];
    os << " " << eVal;
  }
  os << "\n";
}

Face ReadFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  size_t len;
  int eVal;
  is >> len;
  Face eFace(len);
  for (size_t i=0; i<len; i++) {
    is >> eVal;
    eFace[i]=eVal;
  }
  return eFace;
}


vectface ReadListFace(std::istream & is)
{
  if (!is.good()) {
    std::cerr << "ReadListFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  size_t nbFace;
  is >> nbFace;
  if (nbFace == 0) {
    std::cerr << "We cannot handle that case because we need the base length\n";
    throw TerminalException{1};
  }
  Face f = ReadFace(is);
  vectface ListFace(f.size());
  ListFace.push_back(f);
  for (size_t iFace=1; iFace<nbFace; iFace++) {
    Face f2 = ReadFace(is);
    ListFace.push_back(f2);
  }
  return ListFace;
}

void WriteListFace(std::ostream & os, vectface const& ListFace)
{
  size_t nbFace=ListFace.size();
  os << nbFace << "\n";
  for (size_t iFace=0; iFace<nbFace; iFace++)
    WriteFace(os, ListFace[iFace]);
}


void WriteFaceGAP(std::ostream &os, Face const& f)
{
  size_t nb=f.count();
  //  int siz=f.size();
  os << "[";
  boost::dynamic_bitset<>::size_type aPos=f.find_first();
  for (size_t i=0; i<nb; i++) {
    if (i>0)
      os << ",";
    int eVal=int(aPos) + 1;
    os << eVal;
    aPos=f.find_next(aPos);
  }
  os << "]";
}


void WriteListFaceGAP(std::ostream & os, vectface const& ListFace)
{
  os << "[";
  bool IsFirst=true;
  for (auto & eFace : ListFace) {
    if (!IsFirst)
      os << ",";
    IsFirst=false;
    WriteFaceGAP(os, eFace);
  }
  os << "]";
}

void WriteListFaceGAPfile(std::string const& eFile, vectface const& ListFace)
{
  std::ofstream os(eFile);
  os << "return ";
  WriteListFaceGAP(os, ListFace);
  os << ";\n";
}




// We require x and y to be of the same size
bool operator<(Face const& x, Face const& y)
{
  size_t len=x.size();
  for (size_t i=0; i<len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}




void PrintVectInt(std::ostream &os, Face const& eList)
{
  size_t len=eList.size();
  for (size_t i=0; i<len; i++)
    if (eList[i] == 1)
      os << " " << i;
  os << "\n";
}







Face FullFace(size_t const& len)
{
  Face eFace(len);
  for (size_t u=0; u<len; u++)
    eFace[u]=1;
  return eFace;
}



ulong FaceToUnsignedLong(Face const& f)
{
  size_t len=f.size();
  if (len > 32) {
    std::cerr << "Too large value, conversion impossible";
    throw TerminalException{1};
  }
  ulong pos=0;
  ulong pow=1;
  for (size_t i=0; i<len; i++) {
    pos += pow*f[i];
    pow *= 2;
  }
  return pos;
}

Face UnsignedLongToFace(size_t const& len, ulong const& eVal)
{
  if (len > 32) {
    std::cerr << "length error\n";
    throw TerminalException{1};
  }
  ulong eWork = eVal;
  Face eFace(len);
  ulong pow=1;
  for (size_t i=0; i<len; i++) {
    ulong Pow2=pow * 2;
    ulong res=eWork % Pow2;
    if (res == pow) {
      eFace[i]=1;
      eWork -= pow;
    }
    pow=Pow2;
  }
  return eFace;
}

void VectVectInt_Magma_Print(std::ostream &os, vectface const&ListOrbit)
{
  size_t nbOrbit=ListOrbit.size();
  os << "[";
  for (size_t iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    if (iOrbit>0)
      os << ",\n";
    Face eRepr=ListOrbit[iOrbit];
    size_t siz=eRepr.count();
    os << "[";
    boost::dynamic_bitset<>::size_type eVal=eRepr.find_first();
    for (size_t i=0; i<siz; i++) {
      if (i>0)
	os << ",";
      os << int(eVal);
      eVal = eRepr.find_next(eVal);
    }
    os << "]";
  }
  os << "]\n";
}



void VectVectInt_Gap_Print(std::ostream &os, vectface const&ListOrbit)
{
  size_t nbOrbit=ListOrbit.size();
  os << "[";
  for (size_t iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    if (iOrbit>0)
      os << ",\n";
    Face eRepr=ListOrbit[iOrbit];
    size_t siz=eRepr.count();
    os << "[";
    boost::dynamic_bitset<>::size_type eVal=eRepr.find_first();
    for (size_t i=0; i<siz; i++) {
      if (i>0)
	os << ",";
      os << int(eVal+1);
      eVal = eRepr.find_next(eVal);
    }
    os << "]";
  }
  os << "]";
}





#endif
