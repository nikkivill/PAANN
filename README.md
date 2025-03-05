# PAANN

## Description: 
PAANN is a web application for browsing single nucleotide polymorphisms (SNPs) associated with Type-2 Diabetes (T2D) and related phenotypes. 
It provides positive selection statistics for each SNP for Bengali, from Bangladesh and Punjabi from Lahore, Pakistan populations.

### Features:
- Retrieve SNP information by rsID, genomic positions or mapped gene
- Allows for filtering by population or phenotype
- Provides SNP association data information (reference and alternate alleles, mapped gene(s) (if any), phenotype, p-value, ancestry)
- Mapped genes are linked to external websites to provide functional information (PAN-GO/Ensembl GRCh37) 
- Provides population-specific summary statistics of positive selection for SNPs 
- Generates interactive graph visualizations of the summary statistics
- Download data as a .txt file with SNP information and summary statistics values
- Download summary statistics graphs as .png files
  
#### Installation:
1. Clone the git repository: git clone https://github.com/nikkivill/PAANN.git
2. Install dependencies (found as requirements.text): pip install -r requirements.txt
3. Download the dump.sql: git directory PAANN/database 
4. Connect to MySQL command-line: e.g. myysql -u username -p
5. Initialise an empty schema and exit the SQL command-line: e.g. CREATE DATABASE PAANN; exit;
6. Insert dump.sql into empty schema: e.g. mysqldump -u username -p PAANN < dump.sql
7. Run the app.py: python app.py
8. Enter the engine credentials: username, password, host and database name

#### Usage:
- Access the app via the link e.g. '_http://localhost:5000_'
- Enter a search value and optionally filter for population/phenotype
- View results alongside corresponding visualisation of statistics
- Select SNP of interest for further information

#### Programmes used:
- *Database*: MySQL
- *Backend*: FLASK, Python
- *Frontend*: HTML, CSS, JAVASCRIPT
