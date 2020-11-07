TestSingleMatrix:=function(eMat)
  local FileTMP, output, TheCommand, eHNF1, eHNF2;
  FileTMP:="./tmp_in.txt";
  RemoveFileIfExist(FileTMP);
  output:=OutputTextFile(FileTMP, true);
  CPP_WriteMatrix(output, eMat);
  CloseStream(output);
  #
  TheCommand:=Concatenation("./ColHermiteNormalForm ", FileTMP, " out");
  Exec(TheCommand);
  #
  eHNF1:=ReadAsFunction("out")()[2];
  #
  eHNF2:=TransposedMat(HermiteNormalFormIntegerMatTransform(TransposedMat(eMat)).normal);
  if eHNF1 <> eHNF2 then
    Print("eHNF1=", eHNF1, "\n");
    Print("eHNF2=", eHNF2, "\n");
    Print("Inconsistency in the HNF\n");
    Error("Please debug");
  fi;
end;

RandomTestHNF:=function(nbIter, nbRow, nbCol)
  local iIter, eMat, iRow, iCol;
  for iIter in [1..nbIter]
  do
    Print("iIter=", iIter, " / ", nbIter, " nbRow=", nbRow, " nbCol=", nbCol, "\n");
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

SequenceTest:=function()
  RandomTestHNF(100, 4, 3);
  RandomTestHNF(100, 5, 5);
  RandomTestHNF(100, 5, 7);
  RandomTestHNF(100, 7, 5);
  RandomTestHNF(100, 6, 6);
  RandomTestHNF(100, 8, 6);
  RandomTestHNF(100, 6, 8);
  RandomTestHNF(100, 7, 7);
  RandomTestHNF(100, 7, 9);
  RandomTestHNF(100, 10, 10);
  RandomTestHNF(100, 20, 10);
  RandomTestHNF(100, 10, 20);
end;

SequenceTest();
