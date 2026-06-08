#!/bin/bash

#A: Calculate Fst by bcftools.

vcftools --vcf input_data.vcf \
	--weir-fst-pop population_1.txt \
	--weir-fst-pop population_2.txt \
	--out pop1_vs_pop2
#--weir-fst-pop: Specifies the population files for comparison.
#--out: Defines the prefix for output files.

#B: Calculate PBS by Fst.
#Defining a three population tree with focal population A, a relatively closely related population B, and a ‘population outgroup’ C, we can estimate the corresponding branch lengths of a,b,c. Because the pairwise distances T are approximately additive, the branch length a from the three pairwise T values is estimated as:

#T = -log(1 - Fst)

#PBS.A = (Tab + Tac - Tbc)/2

for i in {1..21}; do
    echo "Processing chromosome $i"
    vcftools --gzvcf path/to/chr${i}.vcf.gz \
             --weir-fst-pop sampleA.list \
             --weir-fst-pop sampleB.list \
             --out chr${i}.AvsB
done
