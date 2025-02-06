-- Example Query to list general information for SNPs
SELECT
	sc.SNP_ID,
    sc.Genomic_Coordinates,
    sc.P_Value,
    sc.Mapped_Gene,
    ps.Population_ID,
    a.Reference_Allele,
    a.Alternate_Allele
FROM SNP_central sc
LEFT JOIN Population_SNP ps ON sc.SNP_ID = ps.SNP_ID
LEFT JOIN Population p ON ps.Population_ID = p.Population_ID
LEFT JOIN Allele a ON sc.SNP_ID = a.SNP_ID;

-- Query to find 'Population Specifc Allele Frequency for given SNP'
SELECT
	sc.SNP_ID,
    ps.Population_ID,
    a.Reference_Allele,
    a.Alternate_Allele,
    a.Allele_Frequency
FROM SNP_central sc
LEFT JOIN Population_SNP ps ON sc.SNP_ID = ps.SNP_ID
LEFT JOIN Allele a ON sc.SNP_ID = a.SNP_ID;

-- Query to get 'Summary Statistics of SNPs'
SELECT
	sc.SNP_ID,
    sc.P_Value,
    ps.Population_ID,
	a.Reference_Allele,
    a.Alternate_Allele,
    a.Allele_Frequency,
    ss.Tajima_D,
    ss.FST,
    ss.iHS,
    ss.MK
FROM SNP_central sc
LEFT JOIN Population_SNP ps ON sc.SNP_ID = ps.SNP_ID
LEFT JOIN Allele a ON sc.SNP_ID = a.SNP_ID
LEFT JOIN Summary_Statistics ss ON sc.SNP_ID = ss.SNP_ID

-- Query to 'Find SNPs with significant P-Value'
SELECT
	sc.SNP_ID,
    sc.Genomic_Coordinates,
    sc.P_Value,
    sc.Mapped_Gene,
    ps.Population_ID
FROM SNP_central sc
LEFT JOIN Population_SNP ps ON sc.SNP_ID = ps.SNP_ID
WHERE sc.P_Value < '0.05';

	