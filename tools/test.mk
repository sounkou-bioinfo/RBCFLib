# Compile bcf_field_indexer
bcf_field_indexer: bcf_field_indexer.c
	gcc -O2 -Wall -std=c11 bcf_field_indexer.c -I../src/bcftools-1.22/htslib-1.22 -I. ../src/bcftools-1.22/htslib-1.22/libhts.a -ldeflate -lm -lz -lbz2 -llzma -lcurl -lcrypto -lssl -lpthread -o bcf_field_indexer

simple: test_mmap.c
	gcc -O2 test_mmap.c ../src/hfile_mmap.c -I/usr/local/lib/R/site-library/RBCFLib/include/htslib -I../src/bcftools-1.22/htslib-1.22 -I. ../src/bcftools-1.22/htslib-1.22/libhts.a -ldeflate -lm -lz -lbz2 -llzma -lcurl -lcrypto -lssl -lpthread -o test_mmap

all: bcf_field_indexer test_mmap