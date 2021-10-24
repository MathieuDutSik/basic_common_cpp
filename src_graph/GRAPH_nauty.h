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
  TheCommand=FileConvert + " " + eFileOut.string() + " " + std::to_string(nbVert) + " > " + eFileGrp.string();
  //
  std::ifstream is(eFileGrp.string());
  TheGroupFormat GRP=ReadGroup(is);
  //
  return GRP;
}




template<typename Tgr>
std::option<permlib::Permutation> GRAPH_Isomorphism_Nauty(Tgr const& eGR1, Tgr const& eGR2)
{
  std::string ePrefix=random_string(20);
  TempFile eFileIn ("/tmp/NAUTY_ISOM_" + ePrefix + ".inp");
  TempFile eFileOut("/tmp/NAUTY_ISOM_" + ePrefix + ".out");
  TempFile eFileErr("/tmp/NAUTY_ISOM_" + ePrefix + ".err");
  TempFile eFileIso("/tmp/NAUTY_ISOM_" + ePrefix + ".iso");
  //
  if (eGR1.GetNbVert() != eGR2.GetNbVert() || GetNbColor(eGR1) != GetNbColor(eGR2))
    return {};
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
  return ePerm;
}




#endif
