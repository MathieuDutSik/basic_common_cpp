PROGRAM		= EnumerationCliques_myrvold EnumerationCliques_ordered_enum EnumerationCliques_special_subcliques Block_CppIterator Timing_HashMap Vector_Face_iterator SerializeVectface SerializeFace
#PROGRAM		= EnumerationCliques_myrvold EnumerationCliques_ordered_enum EnumerationCliques_special_subcliques
#PROGRAM		= EnumerationCliques_special_subcliques
#PROGRAM		= Block_CppIterator
#PROGRAM		= Timing_HashMap
#PROGRAM		= Vector_Face_iterator SerializeVectface VectorFace_move
#PROGRAM		= Test_SetCppIterator
#PROGRAM		= Vector_Face_iterator
#PROGRAM		= BlockMultiple_CppIterator
#PROGRAM		= SerializeVectface
#PROGRAM		= SerializeFace
#PROGRAM		= IterFace

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Wall -std=c++17
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= g++ -std=c++17 -Wall -g -O9 -ftemplate-backtrace-limit=0
#CC		= g++ -std=c++17 -Wall -g -ftemplate-backtrace-limit=0 -fsanitize=address -fsanitize=leak -fsanitize=undefined -fbounds-check
# To set in CC_OPTION : -static-libgcc -static-libstdc++
CC		= g++ -std=c++17 -Wall -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CHASHMAP        = -I../sparse-map/include/tsl/ -I../robin-map/include/tsl/ -I../hopscotch-map/include/tsl/
CLOCAL          = -I../src_basic -I../src_number -I../src_matrix -I../src_comb -I../src_xml -I.
CFLAGS 		= -I$(BOOST_INCDIR) -I$(EIGEN_PATH)

CTOTAL 		= $(CLOCAL) $(CHASHMAP) $(CFLAGS)

LIBS		= -lm -pthread $(GMP_CXX_LINK) $(BOOST_LINK)


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
