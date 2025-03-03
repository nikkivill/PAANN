SELECT 
    s.SNP_ID,
    s.RSID,
    s.Chromosome,
    s.Position,
    s.Reference_Allele,
    s.Alternate_Allele,
    CONCAT(
        ROUND(s.P_value / POW(10, FLOOR(LOG10(si.P_value))), 2),'e', FLOOR(LOG10(s.P_value))
	) AS P_value,
    s.Mapped_genes,  
    s.Phenotype,
    COALESCE(st.Population_ID, 'N/A') AS Population,
    cs.FST_Score,
    cs.XPEHH,
    st.Tajima_D,
    st.Haplotype_Diversity,
    st.Nucleotide_Diversity,
    s.Ancestry
FROM 
    SNP_Information s
LEFT JOIN 
    Statistics st ON s.SNP_ID = st.SNP_ID
LEFT JOIN
	Population p ON st.Population_ID = p.Population_ID
LEFT JOIN
	Comparative_Statistics cs ON s.SNP_ID = cs.SNP_ID
