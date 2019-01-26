.PHONY=all

CXX?=g++
CC?=gcc

INCDIRS=. htslib bonsai/clhash/include bonsai bonsai/hll/ bonsai/libpopcnt bonsai/hll/vec bonsai/circularqueue bonsai/pdqsort
INCLUDE=$(patsubst %,-I%,$(INCDIRS))
LD=
LIB=-lcurl -lz

OBJ=bonsai/clhash/clhash.o htslib/libhts.a
all: 10xdash libhts.a
htslib/configure:
	cd htslib && autoreconf
htslib/Makefile: htslib/configure
	cd htslib && ./configure --disable-lzma --disable-bz2
htslib/libhts.a: htslib/Makefile
	+cd htslib && make libhts.a && cp ..
libhts.a: htslib/libhts.a
	cp $< $@
bonsai/clhash/clhash.o: bonsai/clhash/src/clhash.c
	cd bonsai/clhash && make clhash.o
OBJ=libhts.a bonsai/clhash/clhash.o


CXXFLAGS+= -march=native -O3 -std=c++17
FLAGS+= $(CXXFLAGS) -fopenmp

%: src/%.cpp $(OBJ)
	$(CXX)  $(OBJ) $< -o $@ $(INCLUDE) $(FLAGS) $(LD) $(LIB)

