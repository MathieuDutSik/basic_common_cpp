

FCT:=PrintGraphCppPolyhedral_TrivialColor;

PL1:=ArchimedeanPolyhedra("Dodecahedron");
GRA1:=PlanGraphToGRAPE(PL1);
File1:="Dodecahedron.graph";
FCT(File1, GRA1);


PL2:=GoldbergConstruction(PL1, 1, 1);
GRA2:=PlanGraphToGRAPE(PL2);
File2:="Dodecahedron_11.graph";
FCT(File2, GRA2);


PL3:=GoldbergConstruction(PL1, 2, 1);
GRA3:=PlanGraphToGRAPE(PL3);
File3:="Dodecahedron_21.graph";
FCT(File3, GRA3);


PL4:=ArchimedeanPolyhedra("Cube");
GRA4:=PlanGraphToGRAPE(PL4);
File4:="Cube.graph";
FCT(File4, GRA4);


PL4:=ArchimedeanPolyhedra("Tetrahedron");
GRA4:=PlanGraphToGRAPE(PL4);
File4:="Tetrahedron.graph";
FCT(File4, GRA4);
