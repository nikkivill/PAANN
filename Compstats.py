#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# ## Hudson's FST 

# In[1]:
# First install the required packages
pip install scikit-allel pandas numpy pysam

# In[2]:
# Load the packages
import allel
import pandas as pd
import numpy as np

# Then load the files 
beb_vcf = "final_BEB_all_filtered_snps.vcf.gz"
pjl_vcf = "final_PJL_all_filtered_snps.vcf.gz"

# and read them using read_vcf
beb = allel.read_vcf(beb_vcf, fields=['variants/ID', 'variants/POS', 'variants/CHROM', 'calldata/GT'])
pjl = allel.read_vcf(pjl_vcf, fields=['variants/ID', 'variants/POS', 'variants/CHROM', 'calldata/GT'])

# I created a Dataframe of both populations containing ID, chromosome, and position
beb_df = pd.DataFrame({'rsID': beb['variants/ID'], 'CHROM': beb['variants/CHROM'], 'POS': beb['variants/POS']})
pjl_df = pd.DataFrame({'rsID': pjl['variants/ID'], 'CHROM': pjl['variants/CHROM'], 'POS': pjl['variants/POS']})

# merge on common SNPs
common_snps = beb_df.merge(pjl_df, on=['rsID', 'CHROM', 'POS'])

# extract the genotypes
beb_gt = allel.GenotypeArray(beb['calldata/GT'])
pjl_gt = allel.GenotypeArray(pjl['calldata/GT'])

beb_ac = beb_gt.count_alleles()
pjl_ac = pjl_gt.count_alleles()

# compute Hudson's Fst using this built-in formula
num, den = allel.hudson_fst(beb_ac, pjl_ac)

# fix NaN division error so that division by 0 isn't calculated 
fst_values = np.divide(num, den, out=np.zeros_like(num, dtype=float), where=den > 0) 

# put the results into a dataframe
results = pd.DataFrame({
    'rsID': common_snps['rsID'],
    'CHROM': common_snps['CHROM'],
    'POS': common_snps['POS'],
    'Fst': fst_values
})

# then I created DataFrame where each SNP (rsID) from common_snps is associated with two populations ('BEB' and 'PJL') using concat, 
# This basically duplicates each rsID while allocating distinct population labels.
snp_population = pd.concat([
    pd.DataFrame({'rsID': common_snps['rsID'], 'Population': 'BEB'}),
    pd.DataFrame({'rsID': common_snps['rsID'], 'Population': 'PJL'})
], ignore_index=True)

# save the results to csv 
results.to_csv("fst_results.csv", index=False)
snp_population.to_csv("snp_population_table.csv", index=False)

# In[3]:
# Print summary statistics while handling NANs 
print("Median Fst:", np.nanmedian(fst_values))
print("Max Fst:", np.nanmax(fst_values))
print("Min Fst:", np.nanmin(fst_values))

# load the packages for plotting 
import seaborn as sns
import matplotlib.pyplot as plt

# Plot a histogram of Fst values
sns.histplot(fst_values, bins=50, kde=True, color='blue')
plt.xlabel("Fst")
plt.ylabel("Frequency")
plt.title("Fst Distribution Between BEB and PJL")
plt.show()
print("Mean Fst:", np.nanmean(fst_values))

# In[4]:
# load the CSV file
df = pd.read_csv('fst_results.csv')

# round the 'Fst' column to 3 decimal places
df['Fst'] = df['Fst'].round(3)

# Save the modified Dataframe back to a new CSV file
df.to_csv('rounded_fst.csv', index=False)

# In[5]:
# load the CSV file
df = pd.read_csv('rounded_fst.csv')

# rename the columns
df.rename(columns={
    'rsID': 'RSID',
    'CHROM': 'Chromosome',
    'POS': 'Position',
    'Fst': 'FST_Score'
}, inplace=True)

# save the updated CSV file
df.to_csv('fst.csv', index=False)

# ## XPEHH - Cross Population Extended Haplotype Homozygosity 

# In[3]:
pip install scikit-allel numpy pandas

# In[41]:
import allel
from allel.stats.selection import ehh_decay
from scipy.stats import norm
import numpy as np
import pandas as pd

# Firstly load the VCF files
beb_callset = allel.read_vcf('final_BEB_all_filtered_snps.vcf.gz')
pjl_callset = allel.read_vcf('final_PJL_all_filtered_snps.vcf.gz')

# In[42]:
# then I extracted genotype data for the populations from the callset
beb_genotypes = beb_callset['calldata/GT']
pjl_genotypes = pjl_callset['calldata/GT']

# In[43]:
# then I extracted the haplotypes and genotypes
beb_haplotypes = allel.GenotypeArray(beb_genotypes).to_haplotypes()
pjl_haplotypes = allel.GenotypeArray(pjl_genotypes).to_haplotypes()

# and the positions
positions = beb_callset['variants/POS']
print(beb_haplotypes)

# In[46]:
# now to compute the EHH deacy for SNPs, I created two empty lists 
beb_ehh_values = []
pjl_ehh_values = []

# this code loops through the list of SNPs and calculates EHH values for both the haplotypes 
for snp_index in positions:
    beb_ehh = ehh_decay(beb_haplotypes, snp_index)
    pjl_ehh = ehh_decay(pjl_haplotypes, snp_index)
    
    # and I used if-else statement to calculate the mean EHH value is an array or list
    if isinstance(beb_ehh, (list, np.ndarray)):
        beb_ehh_values.append(np.nanmean(beb_ehh))  # mean for BEB
    else:
        beb_ehh_values.append(beb_ehh)  

    if isinstance(pjl_ehh, (list, np.ndarray)):
        pjl_ehh_values.append(np.nanmean(pjl_ehh))  # mean for PJL
    else:
        pjl_ehh_values.append(pjl_ehh) 
        
# then I converted the lists to NumPy arrays
beb_ehh_values = np.array(beb_ehh_values)
pjl_ehh_values = np.array(pjl_ehh_values)

# In[47]:
# avoided division by zero by placing small EHH values with min threshold 
beb_ehh_values[beb_ehh_values < 1e-6] = 1e-6
pjl_ehh_values[pjl_ehh_values < 1e-6] = 1e-6

# then I computed XP-EHH using log-ratio of EHH decay
xpehh_scores = np.log(beb_ehh_values / pjl_ehh_values)

# Normalized the scores using mean and sd
xpehh_score = (xpehh_scores - np.mean(xpehh_scores)) / np.std(xpehh_scores)

# In[48]:
# extract rsID, Chromosome, and Position from data file
rsIDs = beb_callset['variants/ID']  
chromosomes = beb_callset['variants/CHROM']  
positions = beb_callset['variants/POS']  

# finally make a dataframe and save the results 
results = pd.DataFrame({
    'rsID': rsIDs,
    'Chromosome': chromosomes,
    'Position': positions,
    'XPEHH_Normalized': xpehh_scores_normalized
})

results.to_csv('xpehh_scores_final1.csv', index=False)
print(results.head())

# In[50]:
# now load the file 
df = pd.read_csv('xpehh_scores_final1.csv')

# round the 'Fst' column to 3 decimal places
df['XPEHH_Normalized'] = df['XPEHH_Normalized'].round(3)

# save the new DataFrame back to a new CSV file
df.to_csv('XPEHH.csv', index=False)

# In[ ]:

