PROGRAM		= AnalysisSingMat NullspaceComputation ZbasisComputation RowHermiteNormalForm SerializeMatrix TestMatrixInverse TestPerformanceHNF ColHermiteNormalForm NullspaceIntTrMat NullspaceIntMat SolutionIntMat IntersectionLattice TestHilbertMatrix TestProdMat Matrix_As_Set TestNullspaceComputation ScalarCanonicalization MatrixSubsetSolver DeterminantMat ComputePointer NonUniqueScaleToInteger NullspaceComputationExperiment SmithNormalForm TestSubspaceCompletion RemoveFractionMatrix SolutionMatRepetitive TestFindIsotropicMod TestEquaSpace RankMat
#PROGRAM		= TestNullspaceComputation

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= g++ -std=c++17 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -g -O3
#CC		= clang++ -std=c++17 -Wall -Wconversion -g -O3
#CC		= clang++ -std=c++17 -Weverything -g -O3
#CC		= g++ -std=c++17 -Wall -g -O3 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fbounds-check
#CC		= g++ -std=c++17 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_comb


CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) -I.

LIBS		= -lm $(GMP_CXX_LINK) $(BOOST_LINK) -pthread

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
