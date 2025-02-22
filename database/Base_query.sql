SELECT 
    SNP_Information.RSID,
    SNP_Information.Chromosome,
    SNP_Information.Position,
    SNP_Information.Reference_Allele,
    SNP_Information.Alternate_Allele,
    SNP_Information.Ancestry,
	GROUP_CONCAT(Gene.Gene_Name) AS Genes, -- groups all the genes together instead of duplicating entire row 
    CONCAT(
        ROUND(SNP_Information.P_value / POW(10, FLOOR(LOG10(SNP_Information.P_value))), 2),
        'e',
        FLOOR(LOG10(SNP_Information.P_value))
    ) AS P_value,  -- changes the notation of the p-value so that are written in scientific notation with 2 decimal places
    Statistics.Population_ID,
    Statistics.Tajima_D,
    Statistics.Haplotype_Diversity,
    Statistics.Nucleotide_Diversity
FROM 
    SNP_Information
JOIN 
    Statistics ON SNP_Information.SNP_ID = Statistics.SNP_ID
JOIN 
    SNP_Gene ON SNP_Information.SNP_ID = SNP_Gene.SNP_ID 
JOIN 
    Gene ON SNP_Gene.Gene_ID = Gene.Gene_ID  
GROUP BY -- groups the genes based on the information in the same row by aggregating them
    SNP_Information.SNP_ID, 
    SNP_Information.RSID,
    SNP_Information.Chromosome,
    SNP_Information.Position,
    SNP_Information.Reference_Allele,
    SNP_Information.Alternate_Allele,
    SNP_Information.Ancestry,
    Statistics.Population_ID,
    Statistics.Tajima_D,
    Statistics.Haplotype_Diversity,
    Statistics.Nucleotide_Diversity,
    SNP_Information.P_value;

    
