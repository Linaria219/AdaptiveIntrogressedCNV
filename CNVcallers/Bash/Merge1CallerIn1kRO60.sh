#!/bin/bash
ref="path/to/reference.fa"
filelist=`cat software.fileline.txt`
taglist=`cat software.tagsline.txt`
outdir="path/to/output"
sup="path/to/supplement.vcf"  # A VCF file containing the header.

path/to/mergeSVcallers -a ${ref} \
	-f ${filelist} \
	-t ${taglist} \
	-s 1000 \
	-r 0.6 2> ${outdir}/mergedSoftware.err > ${outdir}/mergedSoftware.vcf

awk -v insert="$sup" '
BEGIN {
    inserted = 0
}
/^#CHROM/ && !inserted {
    while ((getline line < insert) > 0)
        print line
    close(insert)
    inserted = 1
}
{ print }
' ${outdir}/mergedSoftware.vcf > ${outdir}/mergedSoftware.sup.vcf

bcftools sort -o ${outdir}/mergedSoftware.srt.vcf ${outdir}/mergedSoftware.sup.vcf
bgzip ${outdir}/mergedSoftware.srt.vcf
tabix -p vcf ${outdir}/mergedSoftware.srt.vcf.gz

