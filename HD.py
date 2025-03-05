## note: match population is ai coded, basically checks file begin name, extracts as pop and puts it in pop column later for sql stuff
import csv
from collections import Counter
import gzip
import re

def calculate_haplotype_frequencies(vcf_file_path):
    match = re.search(r'([A-Z]+)_all_filtered_snps', vcf_file_path)
    population = match.group(1) if match else "Unknown"
    
    output_csv_path = f"{population}_haplotype_frequencies_and_diversity.csv"

    # opens gzip file, writes as csv with the headings wanted
    with gzip.open(vcf_file_path, 'rt') as f, open(output_csv_path, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['SNP_ID', 'Chromosome', 'Position', 'p(0)', 'p(1)', 'Haplotype_Diversity', 'Population'])

        # okay so iterates each line of vcf, uses indexing to find where things are, also skipped # headersss
        for line in f:
            if not line.startswith('##') and not line.startswith('#'):
                columns = line.strip().split('\t')

                # took out some headers so I can put them back into my csv later on, genotype was 9 onwards
                # also initialised haplotypes one so can be added onto CSV later
                chromosome = columns[0]
                position = columns[1]
                snp_id = columns[2] 
                genotypes = columns[9:]
                haplotypes = [] 

                # so splits genotype based on '|', becomes separate alleles which are added to haplotypes list
                # basically creates a haplotype list form initial genotype list
                for genotype in genotypes:
                    allele1, allele2 = genotype.split('|')
                    haplotypes.append(allele1)
                    haplotypes.append(allele2)

                # counts how many unique haplotypes (should be '0' and '1') there are all together in the list
                # Then calcs frequency with formula so basically doing f/2n
                haplotype_counts = Counter(haplotypes)
                total_haplotypes = sum(haplotype_counts.values())
                p_0 = haplotype_counts.get('0', 0) / total_haplotypes
                p_1 = haplotype_counts.get('1', 0) / total_haplotypes
                
                # then we use formula so we N/N-1 * 1-((p_0)^2+(p_1)^2)
                N = total_haplotypes
                diversity = (N / (N - 1)) * (1 - (p_0**2 + p_1**2))

                # add a row on top of the CSV based on what we calc and also previous headers I wanted
                writer.writerow([snp_id, chromosome, position, p_0, p_1, diversity, population])
    
    print(f"Output written to {output_csv_path}")
