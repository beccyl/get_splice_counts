
INC=-I/group/bioi1/nadiad/software/samtools-1.9/ -I/group/bioi1/nadiad/software/samtools-1.9/htslib-1.9
LD=-L/group/bioi1/nadiad/software/samtools-1.9/htslib-1.9 -L/group/bioi1/nadiad/software/samtools/ -L/group/bioi1/nadiad/lib/
HEADERS="BamReader.h"

all: get_splice_counts

get_splice_counts: get_splice_counts.c++ BamReader.h
	g++ -g -std=c++1y -O3 $(INC) -o $@ $< $(HEADERS) $(LD) -lhts -lz -lbz2 -llzma -lcurl #-lprofiler 

clean:                   
	-rm *~ get_splice_counts

