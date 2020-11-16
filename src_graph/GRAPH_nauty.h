#ifndef INCLUDE_GRAPH_NAUTY
#define INCLUDE_GRAPH_NAUTY


template<typename Tgr>
int GetNbColor(Tgr const& eGR)
{
  int nbColor=0;
  int nbVert=eGR.GetNbVert();
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eColor=eGR.GetColor(iVert);
    if (eColor > nbColor)
      nbColor=eColor;
  }
  nbColor++;
  return nbColor;
}


template<typename Tgr>
void NAUTY_PrintPartition(std::ostream & os, Tgr const& eGR)
{
  int nbVert=eGR.GetNbVert();
  int nbColor=GetNbColor(eGR);
  os << "f=[";
  for (int iColor=0; iColor<nbColor; iColor++) {
    bool IsFirst=true;
    for (int iVert=0; iVert<nbVert; iVert++) {
      if (!IsFirst) {
        os << " ";
      }
      IsFirst=false;
      if (eGR.GetColor(iVert) == iColor)
        os << iVert;
    }
    if (iColor < nbColor-1)
      os << "|";
  }
  os << "]\n";
}




template<typename Tgr>
void NAUTY_PrintGraph(std::ostream & os, Tgr const& eGR)
{
  os << "g\n";
  int nbVert=eGR.GetNbVert();
  for (int iVert=0; iVert<nbVert; iVert++) {
    os << iVert << " :";
    for (int jVert=0; jVert<nbVert; jVert++)
      if (eGR.IsAdjacent(iVert, jVert))
        os << " " << jVert;
    os << ";\n";
  }
}



template<typename T, typename Tgr>
void NAUTY_AUTO_WriteFile(std::ostream & os, std::vector<MyMatrix<T>> const& ListMatrix, MyMatrix<T> const& SHV)
{
  WeightMatrix<std::vector<T>, T> ScalMat = GetWeightMatrix_ListMatrix(ListMatrix, SHV);
  Tgr eGR=GetGraphFromWeightedMatrix<std::vector<T>, T, Tgr>(ScalMat);
  int nbVert=eGR.nbVertices;
  os << "n=" << nbVert << "\n";
  NAUTY_PrintPartition(os, eGR);
  NAUTY_PrintGraph(os, eGR);
}



template<typename T, typename Tgr>
bool NAUTY_ISOM_WriteFile(std::ostream & os,
                         std::vector<MyMatrix<T>> const& ListMatrix1, MyMatrix<T> const& SHV1,
                         std::vector<MyMatrix<T>> const& ListMatrix2, MyMatrix<T> const& SHV2)
{
  WeightMatrix<std::vector<T>, T> ScalMat1 = GetWeightMatrix_ListMatrix(ListMatrix1, SHV1);
  WeightMatrix<std::vector<T>, T> ScalMat2 = GetWeightMatrix_ListMatrix(ListMatrix2, SHV2);
  bool test=RenormalizeWeightMatrix(ScalMat1, ScalMat2);
  if (!test)
    return false;
  Tgr eGR1=GetGraphFromWeightedMatrix<std::vector<T>, T, Tgr>(ScalMat1);
  Tgr eGR2=GetGraphFromWeightedMatrix<std::vector<T>, T, Tgr>(ScalMat2);
  int nbVert=eGR1.GetNbVert();
  os << "n=" << nbVert << "\n";
  NAUTY_PrintPartition(os, eGR1);
  NAUTY_PrintGraph(os, eGR1);
  os << "c x @\n";
  NAUTY_PrintGraph(os, eGR2);
  os << "x ##\n";
  return true;
}


template<typename Tgr>
TheGroupFormat GRAPH_Automorphism_Nauty(Tgr const& eGR)
{
  std::string ePrefix=random_string(20);
  TempFile eFileIn ("/tmp/NAUTY_AUTO_" + ePrefix + ".inp");
  TempFile eFileOut("/tmp/NAUTY_AUTO_" + ePrefix + ".out");
  TempFile eFileErr("/tmp/NAUTY_AUTO_" + ePrefix + ".err");
  TempFile eFileGrp("/tmp/NAUTY_AUTO_" + ePrefix + ".grp");
  //
  std::ofstream os(eFileIn.string());
  int nbVert=eGR.GetNbVert();
  os << "n=" << nbVert << "\n";
  NAUTY_PrintPartition(os, eGR);
  NAUTY_PrintGraph(os, eGR);
  os.close();
  //
  std::string FileDR2="dreadnaut";
  std::string TheCommand=FileDR2 + " < " + eFileIn.string() + " > " + eFileOut.string() + " 2> " + eFileErr.string();
  int iret=system(TheCommand.c_str());
  if (iret == -1) {
    std::cerr << "unable to run the command\n";
    throw TerminalException{1};
  }
  if (!eFileOut.IsExisting()) {
    std::cerr << "Error while opening file\n";
    throw TerminalException{1};
  }
  //
  std::string FileConvert="NautyGroupToCPP";
  TheCommand=FileConvert + " " + eFileOut.string() + " " + IntToString(nbVert) + " > " + eFileGrp.string();
  //
  std::ifstream is(eFileGrp.string());
  TheGroupFormat GRP=ReadGroup(is);
  //
  return GRP;
}




template<typename Tgr>
EquivTest<permlib::Permutation> GRAPH_Isomorphism_Nauty(Tgr const& eGR1, Tgr const& eGR2)
{
  std::string ePrefix=random_string(20);
  TempFile eFileIn ("/tmp/NAUTY_ISOM_" + ePrefix + ".inp");
  TempFile eFileOut("/tmp/NAUTY_ISOM_" + ePrefix + ".out");
  TempFile eFileErr("/tmp/NAUTY_ISOM_" + ePrefix + ".err");
  TempFile eFileIso("/tmp/NAUTY_ISOM_" + ePrefix + ".iso");
  //
  if (eGR1.GetNbVert() != eGR2.GetNbVert() || GetNbColor(eGR1) != GetNbColor(eGR2))
    return {false, {}};
  //
  bool test;
  std::ofstream os(eFileIn.string());
  int nbVert=eGR1.GetNbVert();
  os << "n=" << nbVert << "\n";
  NAUTY_PrintPartition(os, eGR1);
  NAUTY_PrintGraph(os, eGR1);
  os << "c x @\n";
  NAUTY_PrintGraph(os, eGR2);
  os << "x ##\n";
  os.close();
  //
  std::string FileDR2="dreadnaut";
  std::string TheCommand=FileDR2 + " < " + eFileIn.string() + " > " + eFileOut.string() + " 2> " + eFileErr.string();
  int iret=system(TheCommand.c_str());
  if (iret == -1) {
    std::cerr << "unable to run the command\n";
    throw TerminalException{1};
  }
  if (!eFileOut.IsExisting()) {
    std::cerr << "Error while opening file\n";
    throw TerminalException{1};
  }
  //
  std::string FileConvert="NautyGroupToCPP";
  TheCommand=FileConvert + " " + eFileOut.string() + " > " + eFileIso.string();
  //
  std::ifstream is(eFileIso.string());
  is >> test;
  if (test == 0) {
    return {false, {}};
  }
  std::vector<permlib::dom_int> eList(nbVert);
  for (int iVert=0; iVert<nbVert; iVert++) {
    int eVal;
    is >> eVal;
    eList[iVert]=eVal;
  }
  permlib::Permutation ePerm(eList);
  return {true, ePerm};
}



template<typename T, typename Tgr>
EquivTest<permlib::Permutation> LinPolytopeGram_Isomorphism_Nauty(std::vector<MyMatrix<T>> const& ListMatrix1, MyMatrix<int> const& SHV1\
, std::vector<MyMatrix<T>> const& ListMatrix2, MyMatrix<int> const& SHV2)
{
  MyMatrix<T> T_SHV1=ConvertMatrixUniversal<T,int>(SHV1);
  MyMatrix<T> T_SHV2=ConvertMatrixUniversal<T,int>(SHV2);
  WeightMatrix<std::vector<T>, T> ScalMat1 = GetWeightMatrix_ListMatrix(ListMatrix1, T_SHV1);
  WeightMatrix<std::vector<T>, T> ScalMat2 = GetWeightMatrix_ListMatrix(ListMatrix2, T_SHV2);
  int nbVert=SHV1.rows();
  bool test=RenormalizeWeightMatrix(ScalMat1, ScalMat2);
  if (!test)
    return {false, {}};
  Tgr eGR1=GetGraphFromWeightedMatrix<std::vector<T>, T, Tgr>(ScalMat1);
  Tgr eGR2=GetGraphFromWeightedMatrix<std::vector<T>, T, Tgr>(ScalMat2);
  EquivTest<permlib::Permutation> eRes=GRAPH_Isomorphism_Nauty(eGR1, eGR2);
  if (!eRes.TheReply)
    return {false, {}};
  std::vector<permlib::dom_int> eListRed(nbVert);
  for (int i=0; i<nbVert; i++)
    eListRed[i]=eRes.TheEquiv.at(i);
  return {true, permlib::Permutation(eListRed)};
}



template<typename T, typename Tgr>
TheGroupFormat LinPolytopeGram_Automorphism_Nauty(std::vector<MyMatrix<T>> const& ListMatrix, MyMatrix<int> const& SHV)
{
  int nbVert=SHV.rows();
  MyMatrix<T> T_SHV=ConvertMatrixUniversal<T,int>(SHV);
  WeightMatrix<std::vector<T>, T> ScalMat = GetWeightMatrix_ListMatrix(ListMatrix, T_SHV);
  Tgr eGR=GetGraphFromWeightedMatrix<std::vector<T>, T, Tgr>(ScalMat);
  TheGroupFormat GRP=GRAPH_Automorphism_Nauty(eGR);
  //
  std::vector<permlib::Permutation> ListGen;
  for (auto & eGen : GRP.group->S) {
    std::vector<permlib::dom_int> v(nbVert);
    for (int iVert=0; iVert<nbVert; iVert++)
      v[iVert]=eGen->at(iVert);
    ListGen.push_back(permlib::Permutation(v));
  }
  return GetPermutationGroup(nbVert, ListGen);
}


#endif
