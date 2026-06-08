#!/bin/bash

sample="YOUR_Sample_Name"
bam_path="path/to/bam"
outdir="path/to/output"
workdir="path/to/workdir"
REF="masked_reference.fa"

ref="${workdir}/hg19_masked.fa"
meta="${workdir}/metadata"

# Step 1: Convert BAM to Fastq and Convert Fastq to 36-mer by using seqkit
samtools view -@ 2 -F 1024 -h "${bam_path}" | samtools fastq -@ 2 - | seqkit sliding -s 36 -W 36 -j 2 | pigz -p 2 > "${outdir}/${sample}.kmer36.fastq.gz"

# Step 2: Use mrsFAST
mrsfast --threads 8 --seq "${outdir}/${sample}.kmer36.fastq.gz" --seqcomp --search "${ref}" -e 2 --disable-nohits -o "${outdir}/${sample}.sam"

# Step 3: Sort BAM
samtools sort -@ 8 -o "${outdir}/${sample}.bam" "${outdir}/${sample}.sam"
samtools index -@ 8 "${outdir}/${sample}.bam"

# Step 4: Calculate depth on standard windows
samtools depth -aa -@ 8 -b "${meta}/standard_windows.gc.bed" "${outdir}/${sample}.bam" > "${outdir}/${sample}_depth_GS.txt"
awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2, $3}' "${outdir}/${sample}_depth_GS.txt" | sort -k1,1 -k2,2n > "${outdir}/${sample}_depth_GS.bed"

# Step 5: GC-depth correction
bedtools map -a "${meta}/standard_windows.gc.bed" -b "${outdir}/${sample}_depth_GS.bed" -c 4 -o mean > "${outdir}/${sample}_standard_gc_depth.bed"
python "${workdir}/bin_gc_depth.py" "${outdir}/${sample}_standard_gc_depth.bed" "${outdir}/${sample}_bin_gc_depth.txt"
pigz -p 22 "${outdir}/${sample}_depth_GS.bed"

# Step 6: CNV per chromosome
CHRS=$(seq 1 22)
export sample
export workdir
export outdir

parallel -j 8 '
chr={};
samtools depth -aa -@ 1 -b ${workdir}/metadata/chr${chr}.windows_unmask.bed ${outdir}/${sample}.bam > ${outdir}/${sample}.chr${chr}.depth.txt && \
awk "BEGIN{{OFS=\"\\t\"}} {{print \$1, \$2-1, \$2, \$3}}" ${outdir}/${sample}.chr${chr}.depth.txt | sort -k1,1 -k2,2n > ${outdir}/${sample}.chr${chr}.depth.bed && \
bedtools map -a ${workdir}/metadata/chr${chr}.windows.bed.gc -b ${outdir}/${sample}.chr${chr}.depth.bed -c 4 -o mean > ${outdir}/${sample}.chr${chr}.gc.depth.bed && \
python ${workdir}/infer_cn.py -d ${outdir}/${sample}_bin_gc_depth.txt -i ${outdir}/${sample}.chr${chr}.gc.depth.bed -o ${outdir}/${sample}.chr${chr}.gc.depth.cn.bed && \
pigz -p 11 ${outdir}/${sample}.chr${chr}.depth.bed && \
' ::: ${CHRS}