#PROGRAM		= GRAPH_EnumerateDistanceSembedding GRAPH_EnumerateShortCycles
#PROGRAM		= GRAPH_CanonicComputation
#PROGRAM		= GRAPH_TestCanonicComputation
#PROGRAM		= GRAPH_TestCanonicComputation GRAPH_CanonicComputation
PROGRAM		= GRAPH_ComputeAutGroup

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++11
#CC		= g++ -std=c++11 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++11 -Wall -g -O2 -static-libgcc -static-libstdc++
CC		= g++ -std=c++11 -Wall -g -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_xml -I.


CFLAGS 		= -I$(GMP_INCDIR)  -I$(PERMLIB_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(MPFR_INCDIR) -I$(EIGEN_PATH) -I$(FLINT_INCDIR) -I$(TBB_INCDIR) -I$(GSL_INCDIR) $(NAUTY_INCLUDE)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(MPFR_LINK) $(FLINT_LINK) $(TBB_LINK) $(GSL_C_LINK) $(NAUTY_LINK)


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
