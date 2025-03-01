#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 8
#$ -l h_rt=08:00:00
#$ -l h_vmem=4G

# load module
module load bcftools

# set input and output directories
INPUT_DIR="../input"
BI_SNPS_DIR="../tmp/bi_snps"
POPULATION_DIR="../tmp/population_snps"

# create output directories
mkdir -p ${BI_SNPS_DIR}
mkdir -p ${POPULATION_DIR}

# loop through chromosome 1-22 VCF files
for INPUT_VCF in $(ls ${INPUT_DIR}/ALL.chr{1..22}*.vcf.gz | sort -V); do    
    # extract chromosome identifier, e.g. "chr1" 
    CHR=$(basename ${INPUT_VCF} | cut -d'.' -f2)

    # set output file name for biallelic SNPs, e.g. "chr1_bi_snps.vcf.gz"
    BI_SNPS_VCF="${BI_SNPS_DIR}/${CHR}_bi_snps.vcf.gz"

    # extract biallelic SNPs and index files for quicker filtering after
    bcftools view --threads 4 --types snps -m2 -M2 -Oz -o ${BI_SNPS_VCF} ${INPUT_VCF}
    bcftools index --threads 4 ${BI_SNPS_VCF}

    # set output file names for populations, e.g. "BEB_chr1_bi_snps.vcf.gz" and "PJL_chr1_bi_snps.vcf.gz"
    BEB_OUTPUT_VCF="${POPULATION_DIR}/BEB_${CHR}_bi_snps.vcf.gz"
    PJL_OUTPUT_VCF="${POPULATION_DIR}/PJL_${CHR}_bi_snps.vcf.gz"

    # filter for Bengali (BEB) population and index files for quicker filtering later
    bcftools view --threads 4 -S BEB_sample_names.txt -Oz -o ${BEB_OUTPUT_VCF} ${BI_SNPS_VCF}
    bcftools index --threads 4 ${BEB_OUTPUT_VCF}

    # filter for Punjabi (PJL) population and index files for quicker filtering later
    bcftools view --threads 4 -S PJL_sample_names.txt -Oz -o ${PJL_OUTPUT_VCF} ${BI_SNPS_VCF}
    bcftools index --threads 4 ${PJL_OUTPUT_VCF}

done




