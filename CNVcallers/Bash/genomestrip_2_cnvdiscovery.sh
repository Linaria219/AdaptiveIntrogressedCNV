#!/bin/bash
#SBATCH --job-name=testGSfor20samples_cnvdiscovery
#SBATCH --partition=YOUR_name_of_partition
#SBATCH --nodes=5
#SBATCH --cpus-per-task=32
#SBATCH -o testGSfor40samples2.out
#SBATCH -e testGSfor40samples2.err

#module load apps/anaconda3/5.2.0
#conda activate genomeSTRiP

which java > /dev/null || exit 1
which Rscript > /dev/null || exit 1
which samtools > /dev/null || exit 1

export SV_DIR="path/to/svtoolkit"
export PATH=${SV_DIR}/bwa:${PATH}
export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH}

classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
rmd="path/to/refmetadata"
rundir="path/to/test_random_40individuals"

java -Xmx16g \
        -cp path/to/svtoolkit/lib/SVToolkit.jar:path/to/svtoolkit/lib/gatk/GenomeAnalysisTK.jar:path/to/svtoolkit/lib/gatk/Queue.jar \
        org.broadinstitute.gatk.queue.QCommandLine \
        -cp path/to/svtoolkit/lib/SVToolkit.jar:path/to/svtoolkit/lib/gatk/GenomeAnalysisTK.jar:path/to/svtoolkit/lib/gatk/Queue.jar \
        -S ${SV_DIR}/qscript/discovery/cnv/CNVDiscoveryPipeline.q \
        -S ${SV_DIR}/qscript/SVQScript.q \
        -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
        -configFile ${SV_DIR}/conf/genstrip_parameters.txt \
        -R ${rmd}/Homo_sapiens_assembly38.fasta \
        -genderMapFile ${rundir}/2183_gender_removeNA.txt \
        -genomeMaskFile ${rmd}/Homo_sapiens_assembly38.mask.101.fasta \
        -ploidyMapFile ${rmd}/Homo_sapiens_assembly38.ploidymap.txt \
	-intervalList ${rmd}/Homo_sapiens_assembly38.interval.list \
        -md ${rundir}/metadata \
        -I ${rundir}/20random.list \
        -tempDir ${rundir}/temp \
        -runDirectory ${rundir} \
        -jobLogDir ${rundir}/cnvdiscovery_joblog \
        -tilingWindowSize 1000 \
        -tilingWindowOverlap 500 \
        -maximumReferenceGapLength 1000 \
        -boundaryPrecision 100 \
        -minimumRefinedLength 500 \
	-jobRunner ParallelShell \
	-gatkJobRunner ParallelShell \
	-maxConcurrentRun 2 \
        -run

