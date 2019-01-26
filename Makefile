CXX?=g++
CC?=gcc

bonsai/clhash/clhash.o: bonsai/clhash/src/clhash.c
	cd bonsai/clhash && make src/clhash.o

INC=-I. -Ihtslib -Ibonsai/clhash/include -Ibonsai -Ibonsai/hll/ -Ibonsai/libpopcnt -Ibonsai/hll/vec -Ibonsai/circularqueue -Ibonsai/pdqsort

%: src/%.cpp htslib/libhts.a bonsai/clhash/clhash.o
	$(CXX)  bonsai/clhash/clhash.o $< -o $@ $(INC) -DUSE_PDQSORT -march=native -lz htslib/libhts.a -lcurl

