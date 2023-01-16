PROGRAM		= TestSquareRoot test_unordered_map_mpzq Factorize PracticalInf TestQuoInt TestTypeConversion RealCubicFieldTest QuadFieldTest RationalTest
#PROGRAM		= test_mpz_givaro
#PROGRAM		= test_unordered_map_mpzq
#PROGRAM		= Factorize
#PROGRAM		= PracticalInf
#PROGRAM		= TestQuoInt
#PROGRAM		= TestConversion
#PROGRAM		= TestSquareRoot
#PROGRAM		= TestTypeConversion
#PROGRAM		= QuadFieldTest
#PROGRAM		= RealCubicFieldTest
#PROGRAM		= RationalTest

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Wall -std=c++17
#CC		= clang++ -g -Wall -std=c++17 -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null
#CC		= clang++ -g -Wall -std=c++17
#CC		= clang++ -g -Wall -std=c++17 -fsanitize=bounds -fsanitize-memory-use-after-dtor -fsanitize-memory-track-origins
#CC		= clang++ -O2 -Wall -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -g -O2 -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
#CC		= g++ -std=c++17 -Wall -O3 -g
#CC		= g++ -std=c++17 -Wall -O3
CC		= g++ -std=c++17 -Wall -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=alignment
#CC		= g++ -std=c++17 -Wall -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o



CLOCAL          = -I../src_basic -I../src_matrix -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(EIGEN_PATH)
CDEBUG 		= -DCHECK_REAL_ALG_NUMERIC
CTOTAL 		= $(CLOCAL) $(CFLAGS) $(CDEBUG)

LIBS		= -lm $(GMP_CXX_LINK)


%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CTOTAL) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
