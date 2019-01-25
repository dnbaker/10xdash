CXX?=g++
CC?=gcc

%: src/%.cpp htslib/libhts.a
	$(CXX) $< -I. -Ihtslib -o $@ -Ibonsai -Ibonsai/hll/ -Ibonsai/libpopcnt -Ibonsai/hll/vec -Ibonsai/circularqueue -Ibonsai/pdqsort -DUSE_PDQSORT -march=native -lz htslib/libhts.a

