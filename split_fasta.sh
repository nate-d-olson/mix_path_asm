#!/usr/bin/sh

for i in {3..17};
	do
		cd NIST_${i}-pilon
		seqretsplit pilon.fasta
		cd ../
	done

