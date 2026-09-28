PROGRAM		= compute_chromatic_number

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++20 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++20
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++20
#CC		= g++ -std=c++20 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++20 -Wall -g -O3
#CC		= g++ -std=c++20 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
CDEPS           = -MMD -MP
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../src_number -I../src_graph -I../src_comb -I../src_matrix


CFLAGS 		= $(MINISAT_INCLUDE) -I$(BOOST_INCDIR) -I$(EIGEN_PATH)  -I$(GMP_INCDIR)

LIBS		= $(MINISAT_CPP_LINK) $(GMP_CXX_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(CDEPS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f *.d
	rm -f $(PROGRAM)

# The header dependencies, written by the compiler as it compiles (-MMD),
# so that a change to a header rebuilds the objects that include it.
-include $(OBJECTS:.o=.d)
