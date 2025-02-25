# PAANN

## Description: 
PAANN is a web application that provides and analyses genomic and functional of **Single Nucleotide Polymorphisms (SNPs)** that are specifically associated with **Type 2 Diabetes (T2D)** associated phenotypes with **positive selection for South Asian populations**.

### Features:
- Retrieve SNP information by rsID, genomic positions or mapped gene(s)
- Provides summary statistics of SNPs in the database
- Generates visualisation plots of the statistics to represent positive selection
- Mapped genes linked to external API to provide KEGG fucntional annotation/ENSEMBL gene page
- Allows for filtering by population or phenotype
- Download data as a text file with SNP information and statistics
  
#### Installation:
1. Clone the git repository: ***{insert link to git repository}***
2. Install dependencies (found as requirments.text): pip install -r requirement.txt
3. Download the dump.sql into SQL workbench: e.g. mysql -u username -p database_name < dump.sql
4. Run the app.py: python app.py
5. Enter the engine credentials: username, password, host and database name

##### Usage:
- Access the app via the link 'http://localhost:5000'
- Enter a search value and optionally filter for population/phenotype
- View results alongside corresponding visualisation of statistics
- Select SNP of interest for further information

###### Technologies/Tools used: EDIT THIS!!!!!!!!!!!!!!!!!!!!!!!!!!
- *Database and data collection*: MySQL, Pandas,
- *Backend*: FLASK, Python
- *Frontend*: HTML, CSS, JAVASCRIPT, DJANGO
- *Visualisation*: Matplotlib, graph.io *****
