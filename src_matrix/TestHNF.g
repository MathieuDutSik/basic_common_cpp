TestSingleMatrix:=function(eMat)
  local FileTMP, output, TheCommand, eHNF1, eHNF2;
  FileTMP:="./tmp_in.txt";
  RemoveFileIfExist(FileTMP);
  output:=OutputTextFile(FileTMP, true);
  CPP_WriteMatrix(output, eMat);
  CloseStream(output);
  #
  TheCommand:=Concatenation("./RowHermiteNormalForm ", FileTMP, " out");
  Exec(TheCommand);
  #
  eHNF1:=ReadAsFunction("out")();
  #
  eHNF2:=HermiteNormalFormIntegerMatTransform(eMat);
  if eHNF1[2]<> eHNF2.normal then
    Print("Inconsistency in the HNF\n");
    Error("Please debug");
  fi;
end;

RandomTestHNF:=function(nbIter, nbRow, nbCol)
  local iIter, eMat, iRow, iCol;
  for iIter in [1..nbIter]
  do
    Print("iIter=", iIter, " / ", nbIter, "\n");
    eMat:=NullMat(nbRow, nbCol);
    for iRow in [1..nbRow]
    do
      for iCol in [1..nbCol]
      do
        eMat[iRow][iCol]:=Random([-3..3]);
      od;
    od;
    TestSingleMatrix(eMat);
  od;
end;
