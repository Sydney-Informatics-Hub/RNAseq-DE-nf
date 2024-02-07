#!/bin/bash
  
#PBS -P er01 
#PBS -N RNA-Seq-nf
#PBS -l walltime=01:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -W umask=022
#PBS -q copyq
#PBS -e RNA-Seq-nf.e
#PBS -o RNA-Seq-nf.o
#PBS -l wd
#PBS -l storage=scratch/er01+gdata/er01
#PBS -l jobfs=4GB

# Load singularity and nextflow modules
module purge
#module load nextflow/22.04.3
module load nextflow/23.04.1  
module load singularity
#module load star/2.7.3a

export NXF_SINGULARITY_CACHEDIR=/scratch/er01/cl9310/singularity

# Fill in these variables for your run
#samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/prep_work/run_pipe/sampleSheet.csv
samples=/scratch/er01/cl9310/1_project/RNAseq-DE-nf/sampleSheet_both_single_paired.csv


#ref=/g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta
#dict=/g/data/er01/SIH-HPC-WGS/Reference/hs38DH.dict
#STARRefIndexPath=/g/data/er01/SIH-Gadi-RNAseq/Reference/GRCh38
#refFasta=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/chromosome22.fasta
#STARRefIndexPath=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/STAR_index_singularity_star_2.7.11a
#refGtf=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/chromosome22.gtf
#SalmonRefIndexPath=
#outDir=results
#NCPUS=2
#adapters_bbmap=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/adapters.fa
#readlen=
#strand="reverse"

# Run the pipeline 
nextflow run main.nf -resume \
        --input ${samples} \
        -profile gadi \
        --adapters_bbmap ${adapters_bbmap}
        --whoami cl9310 --gadi_account er01 \
        --outDir results