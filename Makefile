R_DIR = /usr/local/R-2.15.3
R_INCLUDE = -I$(R_DIR) -I$(R_DIR)/src/include
RMATH_DIR = $(R_DIR)/src/nmath/standalone
RMATH_LINK = $(R_DIR)/src/nmath/standalone/libRmath.a

OPTFLAG = -fprofile-use
PROFFLAG = -fprofile-generate
DEBUGFLAG = -g -O0 
CFLAGS = -Wall -O3

HEADER_FILES = codonString.h

epiPopOpt:
	g++ -o epiPop codonString.cpp seed_time.cpp epiPop.cpp $(CFLAGS) $(OPTFLAG) $(R_INCLUDE) $(RMATH_LINK)

epiPop: codonString.o seed_time.o epiPop.cpp
	g++ -o $@ $^ $(CFLAGS) $(PROFFLAG) $(R_INCLUDE) $(RMATH_LINK)

%.o: %.cpp $(HEADER_FILES)
	g++ -c $(CFLAGS) $(PROFFLAG) $(R_INCLUDE) $(RMATH_LINK) $< 

.PHONY: clean wipe

clean:
	rm -f *.o

wipe:
	rm -f *.o *.gcno *.gcda

