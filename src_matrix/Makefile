PROGRAM		= AnalysisSingMat NullspaceComputation ZbasisComputation RowHermiteNormalForm SerializeMatrix Test_MatrixInverse Test_PerformanceHNF ColHermiteNormalForm NullspaceIntTrMat NullspaceIntMat SolutionIntMat IntersectionLattice Test_HilbertMatrix Test_ProdMat Matrix_As_Set Test_NullspaceComputation ScalarCanonicalization MatrixSubsetSolver DeterminantMat DeterminantMatMod DeterminantMatHadamard ComputePointer NonUniqueScaleToInteger NullspaceComputationExperiment SelectRowColMatModHadamard SmithNormalForm SmithNormalFormIntegerMat Test_SubspaceCompletion RemoveFractionMatrix SolutionMatRepetitive Test_FindIsotropicMod Test_EquaSpace RankMat Matrix_As_SetPosScal NullspaceMatMod IsotropicMatMod RescaleInverse SmithNormalFormInvariant ComputeChromaticLowerBounds

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= clang++ -g -Weverything -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
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
