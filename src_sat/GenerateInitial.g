PL:=ArchimedeanPolyhedra("Dodecahedron");

GRA:=PlanGraphToGRAPE(PL);


GenerateInputFile:=function(GRA)
  local output, i, LAdj, eVal;
  output:=OutputTextFile("InputGraph.txt", true);
  AppendTo(output, GRA.order, "\n");
  for i in [1..GRA.order]
  do
    LAdj:=Adjacency(GRA, i);
    AppendTo(output, Length(LAdj));
    for eVal in LAdj
    do
      AppendTo(output, " ", eVal-1);
    od;
    AppendTo(output, "\n");
  od;
  CloseStream(output);
end;

GenerateInputFile(GRA);
