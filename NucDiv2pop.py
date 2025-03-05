import vcf
import csv

def calculate_nucleotide_diversity(genotypes):
    """
    Calculate nucleotide diversity (π) for a list of genotypes.
    :param genotypes: List of genotypes (strings) for a particular SNP
    :return: Nucleotide diversity (π) value
    """
    n = len(genotypes)
    if n < 2:
        return 0  # Less than 2 genotypes means no diversity

    # Count pairwise differences
    pairwise_diff_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            # If genotypes differ, increment the difference count
            if genotypes[i] != genotypes[j]:
                pairwise_diff_count += 1

    # Nucleotide diversity calculation: average pairwise differences
    return pairwise_diff_count / (n * (n - 1) / 2)

def parse_vcf(file_path):
    """
    Parse a VCF file and return a list of SNPs and their genotypes.
    :param file_path: Path to the VCF file
    :return: List of tuples containing the RSID and list of genotypes
    """
    vcf_reader = vcf.Reader(filename=file_path)
    snp_data = []

    # Iterate through each record in the VCF file
    for record in vcf_reader:
        rsid = record.ID  # SNP ID (e.g., rsID)
        genotypes = []
        
        # For each sample, collect the genotype at this SNP
        for sample in record.samples:
            genotype = sample['GT']
            if genotype is None:
                genotypes.append('NA')  # Missing data (you can skip or handle differently)
            else:
                genotypes.append(''.join(map(str, genotype)))  # Concatenate alleles
        
        snp_data.append((rsid, genotypes))
    
    return snp_data

def calculate_and_output_diversity_multiple(input_vcf_files, output_csv):
    """
    Calculate nucleotide diversity for each SNP across multiple VCF files (populations)
    and output to a single CSV.
    :param input_vcf_files: List of input VCF file paths
    :param output_csv: Path to the output CSV file
    """
    # Initialize a dictionary to store SNP data across all populations
    all_snps = {}

    # Loop through each VCF file and gather SNP data
    for vcf_file in input_vcf_files:
        snp_data = parse_vcf(vcf_file)
        
        # Store the data for each population
        population_name = vcf_file.split('.')[0]  # Extract population name from the VCF file name
        
        for rsid, genotypes in snp_data:
            if rsid not in all_snps:
                all_snps[rsid] = {}
            
            # Store genotypes for each population
            all_snps[rsid][population_name] = genotypes

    # Calculate nucleotide diversity for each SNP
    results = []

    for rsid, populations in all_snps.items():
        diversity_data = {'SNP_ID': rsid}
        for population, genotypes in populations.items():
            # Calculate nucleotide diversity for this SNP in this population
            diversity = calculate_nucleotide_diversity(genotypes)
            diversity_data[population] = diversity
        results.append(diversity_data)

    # Write results to CSV
    with open(output_csv, mode='w', newline='') as file:
        fieldnames = ['SNP_ID'] + list(all_snps[next(iter(all_snps))].keys())  # Use population names as columns
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"Results saved to {output_csv}")

# Example Usage
if __name__ == '__main__':
    input_vcf_files = ['PJL_all_filtered_snps.vcf', 'BEB_all_filtered_snps.vcf']  # List your VCF files here
    output_csv = 'nucleotide_diversity_results.csv'  # Output CSV file path
    calculate_and_output_diversity_multiple(input_vcf_files, output_csv)
