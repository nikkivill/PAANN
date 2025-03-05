import os
import gzip
import sys
import allel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def process_vcf(vcf_file, output_csv, window_size=5000, step_size=2500):
    """ 
    processes VCF file, calculates tajima's D per chromosome, and save the results in a single CSV.
    """
    is_gzipped = vcf_file.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    mode = 'rt' if is_gzipped else 'r'

    print(f"Processing {vcf_file}...")

    callset = allel.read_vcf(vcf_file)

    snp_rsids = callset.get('variants/ID')
    snp_positions = callset['variants/POS']
    snp_chromosomes = callset['variants/CHROM']
    snp_rsids = np.where(snp_rsids == None, "NA", snp_rsids)

    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    allele_counts = genotypes.to_allele_counts()

    print("calculating Tajima's D...")
    tajima_results = []
    
    for chrom in np.unique(snp_chromosomes):
        mask = snp_chromosomes == chrom
        chrom_positions = snp_positions[mask]
        chrom_rsids = snp_rsids[mask]
        chrom_ac = allele_counts.compress(mask, axis=0)

        if chrom_positions.size == 0:
            continue  #skip empty chromosomes

        tajima_values, window_positions, snp_to_tajima = [], [], {}

        for start in range(chrom_positions.min(), chrom_positions.max() - window_size, step_size):
            end = start + window_size
            window_mask = (chrom_positions >= start) & (chrom_positions < end)
            window_ac = chrom_ac.compress(window_mask, axis=0)

            if window_ac.shape[0] > 1:  #ensure at least 2 SNPs
                window_ac = window_ac.sum(axis=1)  #convert from (X, Y, 2) to (X, 2)
                tajima_d = allel.tajima_d(window_ac)

                tajima_values.append(tajima_d)
                window_positions.append((start + end) / 2)

                for pos, rsid in zip(chrom_positions[window_mask], chrom_rsids[window_mask]):
                    snp_to_tajima[rsid] = tajima_d

        df_snp = pd.DataFrame({
            'Chromosome': chrom,
            'RSID': chrom_rsids,
            'Position': chrom_positions,
            'TajimaD': [snp_to_tajima.get(rsid, np.nan) for rsid in chrom_rsids]
        })

        tajima_results.append(df_snp)

    #merge all chromosomes and sort by chromosome and position
    final_df = pd.concat(tajima_results, ignore_index=True)

    #convert chromosome to numeric if possible, otherwise leave as string
    try:
        final_df["Chromosome"] = pd.to_numeric(final_df["Chromosome"])
    except ValueError:
        pass  #keep as string if there are non-numeric chromosomes

    final_df = final_df.sort_values(by=["Chromosome", "Position"], ascending=[True, True])

    #ensure sorting by chromosome  and position
    final_df = final_df.sort_values(by=['Chromosome', 'Position'], key=lambda col: pd.to_numeric(col, errors='coerce'))

    #reorder columns
    final_df = final_df[['RSID', 'Chromosome', 'Position', 'TajimaD']]

    #save CSV
    final_df.to_csv(output_csv, index=False)
    print(f"Final results saved to {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <final_results_csv>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    final_csv_output = sys.argv[2]

    process_vcf(input_vcf, final_csv_output)
