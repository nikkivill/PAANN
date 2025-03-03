# PAANN

## Description: 
PAANN is a web application that identifies, analyses and interprets **Single Nucleotide Polymorphisms (SNPs)** tlinked to **Type 2 Diabetes (T2D)**. Specfically, it focuses on SNPs that exhibit **positive selection for South Asian populations**.

### Features:
- Retrieve SNP information by rsID, genomic positions or mapped gene(s)
- Provides summary statistics of SNPs in the database
- Generates visualisation plots of the statistics to represent positive selection
- Mapped genes linked to external API to provide KEGG fucntional annotation/ENSEMBL gene page
- Allows for filtering by population or phenotype
- Download data as a text file with SNP information and statistics
  
#### Installation:
1. Clone the git repository: https://github.com/nikkivill/PAANN.git
2. Install dependencies (found as requirments.text): pip install -r requirements.txt
3. Download the dump.sql: e.g. mysql -u username -p database_name < dump.sql
4. Create a new schema in MySQL, import the dump.sql and execute the query to populate the database
5. Run the app.py: python app.py
6. Enter the engine credentials: username, password, host and database name

#### Usage:
- Access the app via the link e.g. '_http://localhost:5000_'
- Enter a search value and optionally filter for population/phenotype
- View results alongside corresponding visualisation of statistics
- Select SNP of interest for further information

#### Programmes used:
- *Database*: MySQL
- *Backend*: FLASK, Python
- *Frontend*: HTML, CSS, JAVASCRIPT
