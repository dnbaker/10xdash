.PHONY=all

CXX?=g++
CC?=gcc

INCDIRS=. htslib bonsai/clhash/include bonsai bonsai/hll/ bonsai/libpopcnt bonsai/hll/vec bonsai/circularqueue bonsai/pdqsort
INCLUDE=$(patsubst %,-I%,$(INCDIRS))
LD=
LIB=-lz

all: 10xdash htslib/libhts.a
htslib/libhts.a:
	cd htslib && autoheader && autoconf && \
    ./configure --disable-lzma --disable-bz2 --disable-libcurl && make libhts.a
libhts.a: htslib/libhts.a
	cp htslib/libhts.a libhts.a
bonsai/clhash/clhash.o: bonsai/clhash/src/clhash.c
	cd bonsai/clhash && make clhash.o
OBJ=libhts.a bonsai/clhash/clhash.o


CXXFLAGS+= -march=native -O3 -std=c++14 -fno-strict-aliasing
HEADERS=$(wildcard include/*.h)
WARNINGS=-Wextra -Wall -pedantic -Wno-ignored-attributes -Wno-char-subscripts \
		 -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
		 -Wformat -Wcast-align \
		 -pedantic -Wunused-variable -Wno-attributes -Wno-unused-parameter -Wno-unused-function -Wno-unused-label
        
FLAGS+= $(WARNINGS) $(CXXFLAGS) -fopenmp -DNOT_THREADSAFE -DENABLE_COMPUTED_GOTO -mpclmul -pipe

%: src/%.cpp $(OBJ) $(HEADERS)
	$(CXX)  $< -o $@ $(INCLUDE) $(FLAGS) $(LD) $(LIB) $(OBJ)

clean:
	rm -f 10xdash $(OBJ) htslib/libhts.a
