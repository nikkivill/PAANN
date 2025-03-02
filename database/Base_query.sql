SELECT 
	SNP_Information.RSID,
	SNP_Information.Chromosome,
	SNP_Information.Position,
	SNP_Information.Alternate_Allele,
	SNP_Information.Reference_Allele,
	SNP_Information.Mapped_genes,
	SNP_Information.Phenotype,
	CONCAT(
	ROUND(SNP_Information.P_value / POW(10, FLOOR(LOG10(SNP_Information.P_value))), 2),
	'e',
	FLOOR(LOG10(SNP_Information.P_value))
	) AS P_value,  -- changes the notation of the p-value so that are written in scientific notation with 2 decimal places
	SNP_Information.Ancestry,
	Statistics.Population_ID,
	Comparative_Statistics.Tajima_D,
	Comparative_Statistics.Haplotype_Diversity,
	Comparative_Statistics.Nucleotide_Diversity
FROM 
	SNP_Information
LEFT JOIN 
	Statistics ON SNP_Information.SNP_ID = Statistics.SNP_ID
LEFT JOIN 
	Comparative_Statistics ON SNP_Information.SNP_ID = Comparative_Statistics.SNP_ID
