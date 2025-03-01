#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=02:00:00
#$ -l h_vmem=4G

# load module
module load bcftools

# set input and results directories
input_dir="./population_snps"
results_dir="../final_results"

# temporary directory for filtered chromosome files (will be deleted later after merging)
temp_dir=$(mktemp -d)

# add echos so troubleshooting will be easier 
echo "Starting SNP filtering for each population chromosome..."

# filter for our rsids from association data "all_rsids.txt"
# loop through biallelic snps for BEB and PJL chromosomes 1 to 22
for chr in {1..22}; do
  # filter and save VCFs to temporary directory
  bcftools view --threads 4 -i "ID=@all_rsids.txt" -Oz \
    -o ${temp_dir}/BEB_chr${chr}_filtered.vcf.gz ${input_dir}/BEB_chr${chr}_bi_snps.vcf.gz

  bcftools view --threads 4 -i "ID=@all_rsids.txt" -Oz \
    -o ${temp_dir}/PJL_chr${chr}_filtered.vcf.gz ${input_dir}/PJL_chr${chr}_bi_snps.vcf.gz

  # index the filtered files
  bcftools index --threads 4 ${temp_dir}/BEB_chr${chr}_filtered.vcf.gz
  bcftools index --threads 4 ${temp_dir}/PJL_chr${chr}_filtered.vcf.gz

done

echo "Merging all filtered files into final population VCFs..."

# merge all filtered chromosomes into one final VCF per population
bcftools concat --threads 4 --allow-overlaps -Oz \
  -o ${results_dir}/final_BEB_all_filtered_snps.vcf.gz ${temp_dir}/BEB_chr*_filtered.vcf.gz

bcftools concat --threads 4 --allow-overlaps -Oz \
  -o ${results_dir}/final_PJL_all_filtered_snps.vcf.gz ${temp_dir}/PJL_chr*_filtered.vcf.gz

# index the final merged VCFs
bcftools index --threads 4 ${results_dir}/final_BEB_all_filtered_snps.vcf.gz
bcftools index --threads 4 ${results_dir}/final_PJL_all_filtered_snps.vcf.gz

# remove temporary files
rm -r ${temp_dir}

echo "Filtered and merged files are saved in: ${results_dir}"


