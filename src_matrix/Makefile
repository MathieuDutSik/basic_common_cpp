#PROGRAM		= AnalysisSingMat NullspaceComputation ZbasisComputation RowHermiteNormalForm
#PROGRAM		= ComputePointer
#PROGRAM		= RowHermiteNormalForm
#PROGRAM		= TestMatrixInverse
#PROGRAM		= TestPerformanceHNF
#PROGRAM		= TestQuoInt
PROGRAM		= NullspaceComputationCheck
#PROGRAM		= RowHermiteNormalForm TestPerformanceHNF ColHermiteNormalForm

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
#CC		= clang++ -g -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++11
#CC		= g++ -std=c++11 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++11 -Wall -g -O3
#CC		= g++ -std=c++11 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_comb


CFLAGS 		= -I$(GMP_INCDIR)  -I$(PERMLIB_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH)

LIBS		= -lm $(GMP_CXX_LINK) -pthread

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
