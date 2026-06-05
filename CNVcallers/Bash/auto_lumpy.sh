#!/bin/bash
# Workpath=/public/group_data_2023/caizl/b_CNVcalling/lumpy
current_date=$(date +"%Y%m%d")
filename=$1
filenamebase=$(basename "$filename" .bam)

# Preprocess
# Speedseq align -R "@RG\tID:id\tSM:sample\tLB:lib" \
#       ref.fa \
#       sample.r1.fq \
#       sample.r2.fq

# Extract the discordant paired-end alignments.
samtools view -b -@ 10 -F 1294 ${filename} > discordants.unsorted.${current_date}.${filename}

# Extract the split-read alignments
samtools view -h ${filename} \
	| extractSplitReads_BwaMem -i stdin \
	| samtools view -Sb - \
	> splitters.unsorted.${current_date}.${filename}

samtools sort -@ 20 -o splitters.${current_date}.${filename} splitters.unsorted.${current_date}.${filename}
samtools sort -@ 20 -o discordants.${current_date}.${filename} discordants.unsorted.${current_date}.${filename}
rm splitters.unsorted.${current_date}.${filename}
rm discordants.unsorted.${current_date}.${filename}

# What the hell is: samtools view -r readgroup1 ${filename} \ (in Github)
samtools view -@ 10 ${filename} \
	| tail -n+100000 \
	| pairend_distro.py \
	-r 101 \
	-X 4 \
	-N 10000 \
	-o ${filenamebase}.lib1.histo > statistics_${filenamebase}
samplemean=`awk '{split($1, a, ":"); print int(a[2])}' statistics_${filenamebase}`
samplestdev=`awk '{split($2, a, ":"); print int(a[2])}' statistics_${filenamebase}`

# Run LUMPY (traditional): conda activate lumpy
lumpy \
	-mw 4 \
	-tt 0 \
	-pe id:sample,bam_file:discordants.${current_date}.${filename},histo_file:${filenamebase}.lib1.histo,mean:${samplemean},stdev:${samplestdev},read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
	-sr id:sample,bam_file:splitters.${current_date}.${filename},back_distance:10,weight:1,min_mapping_threshold:20 \
	> ${filenamebase}_${current_date}.vcf

# Call genotypes on LUMPY output VCF files
svtyper \
	-B ${filename} \
	-S splitters.${current_date}.${filename} \
	-i ${filenamebase}_${current_date}.vcf 
	> ${filenamebase}_${current_date}.gt.vcf

