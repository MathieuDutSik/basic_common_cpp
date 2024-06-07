Print("Beginning TestQuadResidue\n");

TestQuadResidue:=function(eRec)
    local a, p, FileOut, eProg, TheCommand, U;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    #
    eProg:="../../src_number/ComputeQuadraticResidue";
    a:=String(eRec[1]);
    p:=String(eRec[2]);
    Print("a=", a, " p=", p, "\n");
    TheCommand:=Concatenation(eProg, " ", a, " ", p, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    if Length(eRec) = 3 then
        if U<>eRec[3] then
            Print("Inconsistency in the values\n");
            return rec(is_correct:=false);
        fi;
    fi;
    return rec(is_correct:=true);
end;

ListRec:=ReadAsFunction("ListQuadTest")();;


FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        RecReply:=TestQuadResidue(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

