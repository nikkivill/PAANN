from flask import Flask, render_template, request, send_file, Response, redirect
import pandas as pd
from sqlalchemy import create_engine, text
import logging
import math
import io
import requests  # Add this import for KEGG API requests

app = Flask(__name__)

# Logging setup
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_db_engine(user, password, host, database):
    """
    Creates and returns a SQLAlchemy engine for MySQL.
    """
    return create_engine(f"mysql+mysqlconnector://{user}:{password}@{host}/{database}")

user = input("Enter MySQL username: ")
password = input("Enter MySQL password: ")
host = input("Enter MySQL host: ")
database = input("Enter MySQL database: ")

engine = create_db_engine(user, password, host, database)


# Test database connection
def test_db_connection():
    if engine is None:
        logger.error("Database engine is not initialized.")
        return False

    try:
        with engine.connect() as conn:
            logger.info("Successfully connected to the database.")
            return True
    except Exception as e:
        logger.error(f"Error connecting to the database: {e}")
        return False


# Home route
@app.route("/")
def home():
    return render_template("home.html")

# Search route
@app.route("/search", methods=["GET", "POST"])
def search():
    if request.method == "POST":
        query = request.form.get("query")
        search_type = request.form.get("filter")
        populations = request.form.getlist("Population")
        phenotypes = request.form.getlist("phenotype")
        chromosome = request.form.get("chromosome")
        start_position = request.form.get("start_position")
        end_position = request.form.get("end_position")
    else:
        query = request.args.get("query")
        search_type = request.args.get("filter")
        populations = request.args.getlist("Population")
        phenotypes = request.args.getlist("phenotype")
        chromosome = request.args.get("chromosome")
        start_position = request.args.get("start_position")
        end_position = request.args.get("end_position")

    page = request.args.get('page', 1, type=int)
    per_page = 20

    if search_type == "genomic_region":
        search_value = None
    else:
        search_value = query + "%" if query else None

    total_results, results = execute_query(search_type, search_value, populations, phenotypes, chromosome, start_position, end_position, page, per_page)
    
    # Combine results with summary stats from SS_PJL, SS_BEB, and Comparative_Statistics
    combined_results = []
    with engine.connect() as connection:
        for result in results:
            rsid = result['RSID']

            # Fetch summary stats from SS_PJL
            summary_pjl = connection.execute(text("""
            SELECT TajimaD, Nucleotide_Diversity, Haplotype_Diversity 
            FROM Pop_Statistics WHERE RSID = :rsid AND Population_ID = 'PJL'
            """), {'rsid': rsid}).fetchone()

            # Fetch summary stats from SS_BEB
            summary_beb = connection.execute(text("""
            SELECT TajimaD, Nucleotide_Diversity, Haplotype_Diversity 
            FROM Pop_Statistics WHERE RSID = :rsid AND Population_ID = 'BEB'
            """), {'rsid': rsid}).fetchone()

            # Fetch comparative statistics from Comparative_Statistics
            comparative_stats = connection.execute(text("""
            SELECT FST_Score, XPEHH_Score 
            FROM Comparative_Statistics WHERE RSID = :rsid
            """), {'rsid': rsid}).fetchone()

            # Convert results to dictionaries if not None
            summary_pjl_dict = dict(summary_pjl._mapping) if summary_pjl else {}
            summary_beb_dict = dict(summary_beb._mapping) if summary_beb else {}
            comparative_dict = dict(comparative_stats._mapping) if comparative_stats else {}

            # Combine result with summary stats
            combined_result = {
                **result,  # Unpack original result
                'TajimaD_PJL': summary_pjl_dict.get('TajimaD'),
                'Nucleotide_Diversity_PJL': summary_pjl_dict.get('Nucleotide_Diversity'),
                'Haplotype_Diversity_PJL': summary_pjl_dict.get('Haplotype_Diversity'),
                'TajimaD_BEB': summary_beb_dict.get('TajimaD'),
                'Nucleotide_Diversity_BEB': summary_beb_dict.get('Nucleotide_Diversity'),
                'Haplotype_Diversity_BEB': summary_beb_dict.get('Haplotype_Diversity'),
                'FST_Score': comparative_dict.get('FST_Score'),
                'XPEHH_Score': comparative_dict.get('XPEHH_Score'),
            }
            combined_results.append(combined_result)

    total_pages = math.ceil(total_results / per_page)
    logger.info(f"Results: {combined_results}")

    return render_template(
        "search_results.html",
        query=query,
        results=combined_results,
        search_type=search_type,
        page=page,
        total_pages=total_pages,
        chromosome=chromosome,
        start_position=start_position,
        end_position=end_position,
        populations=populations,
        phenotypes=phenotypes
    )

def execute_query(search_type, search_value, populations=None, phenotypes=None, chromosome=None, start_position=None, end_position=None, page=1, per_page=20):
    query_map = {
        "snp_id": "s.RSID LIKE :search_value",
        "gene": "g.Mapped_Gene LIKE :search_value",
        "position": "s.Position LIKE :search_value",
        "chromosome": "s.Chromosome LIKE :search_value",
        "genomic_region": "s.Chromosome = :chromosome AND s.Position BETWEEN :start_position AND :end_position"
    }

    if search_type not in query_map:
        logger.error(f"Invalid search type: {search_type}")
        return 0, []

    sql_query = f"""
        SELECT
            s.SNP_ID,
            s.RSID,
            s.Chromosome,
            s.Position,
            s.Reference_Allele,
            s.Alternate_Allele,
            s.P_value,
            g.Mapped_Gene,
            ph.Phenotype AS Phenotype,
            COALESCE(p.Population_ID, 'N/A') AS Population,
            cs.FST_Score,
            cs.XPEHH_Score,
            ps.TajimaD,
            ps.Haplotype_Diversity,
            ps.Nucleotide_Diversity,
            s.Ancestry
        FROM SNP s
        LEFT JOIN Gene g ON s.RSID = g.RSID
        LEFT JOIN Phenotype ph ON s.SNP_ID = ph.SNP_ID
        LEFT JOIN SNP_Population sp ON s.SNP_ID = sp.SNP_ID
        LEFT JOIN Population p ON sp.Population_ID = p.Population_ID
        LEFT JOIN Comparative_Statistics cs ON s.SNP_ID = cs.SNP_ID
        LEFT JOIN Pop_Statistics ps ON s.SNP_ID = ps.SNP_ID AND ps.Population_ID = p.Population_ID
        WHERE {query_map[search_type]}
    """

    params = {"search_value": f"%{search_value}%"} if search_type != "genomic_region" else {
        "chromosome": chromosome,
        "start_position": start_position,
        "end_position": end_position
    }

    if populations:
        placeholders = ",".join([f":pop{i}" for i in range(len(populations))])
        sql_query += f" AND sp.Population_ID IN ({placeholders})"
        for i, pop in enumerate(populations):
            params[f"pop{i}"] = pop


    if phenotypes:
        sql_query += f" AND ph.Phenotype IN ({','.join([':ph' + str(i) for i in range(len(phenotypes))])})"
        for i, ph in enumerate(phenotypes):
            params[f"ph{i}"] = ph

    try:
        with engine.connect() as conn:
            logger.info(f"Executing query: {sql_query}")
            logger.info(f"Query parameters: {params}")

            # Get total results count
            count_query = f"SELECT COUNT(*) FROM ({sql_query}) AS subquery"
            total_results = pd.read_sql(text(count_query), conn, params=params).iloc[0, 0]

            # Pagination
            sql_query += " LIMIT :limit OFFSET :offset"
            params["limit"] = per_page
            params["offset"] = (page - 1) * per_page

            # Fetch results
            results = pd.read_sql(text(sql_query), conn, params=params)
            logger.info(f"Query returned {len(results)} rows.")

            return total_results, results.to_dict(orient='records')
    except Exception as e:
        logger.error(f"Error executing query: {e}")
        return 0, []

    
###MUST GO OVER, THIS IS AI CODED!!!! WASNT WORKING BEFORE COS IT WAS NOT SPECIFIC BUT NOW WORKS, NEED TO UNDERSTAND IT THO
@app.route("/snp/<snp_id>")
def snp_detail(snp_id):
    try:
        with engine.connect() as conn:
            # Fetch basic SNP information
            snp_query = text("""
                SELECT
                    s.SNP_ID,
                    s.RSID,
                    s.Chromosome,
                    s.Position,
                    s.Reference_Allele,
                    s.Alternate_Allele,
                    s.P_value,
                    g.Mapped_Gene,
                    ph.Phenotype AS Phenotype
                FROM SNP s
                LEFT JOIN Gene g ON s.RSID = g.RSID
                LEFT JOIN Phenotype ph ON s.SNP_ID = ph.SNP_ID
                WHERE s.RSID = :snp_id
            """)

            
            snp_data = conn.execute(snp_query, {'snp_id': snp_id}).fetchone()

            if not snp_data:
                return render_template("snp_detail.html", snp=None)

            # Convert the result to a dictionary
            snp = dict(snp_data._mapping)

            # Fetch BEB-specific statistics
            beb_query = text("""
                SELECT
                    TajimaD,
                    Nucleotide_Diversity,
                    Haplotype_Diversity
                FROM Pop_Statistics
                WHERE RSID = :snp_id AND Population_ID = 'BEB'
            """)
            beb_stats = conn.execute(beb_query, {'snp_id': snp_id}).fetchone()
            if beb_stats:
                snp.update({
                    'TajimaD_BEB': beb_stats.TajimaD,
                    'Nucleotide_Diversity_BEB': beb_stats.Nucleotide_Diversity,
                    'Haplotype_Diversity_BEB': beb_stats.Haplotype_Diversity
                })
            else:
                snp.update({
                    'TajimaD_BEB': 'N/A',
                    'Nucleotide_Diversity_BEB': 'N/A',
                    'Haplotype_Diversity_BEB': 'N/A'
                })

            # Fetch PJL-specific statistics
            pjl_query = text("""
                SELECT
                    TajimaD,
                    Nucleotide_Diversity,
                    Haplotype_Diversity
                FROM Pop_Statistics
                WHERE RSID = :snp_id AND Population_ID = 'PJL'
            """)
            pjl_stats = conn.execute(pjl_query, {'snp_id': snp_id}).fetchone()
            if pjl_stats:
                snp.update({
                    'TajimaD_PJL': pjl_stats.TajimaD,
                    'Nucleotide_Diversity_PJL': pjl_stats.Nucleotide_Diversity,
                    'Haplotype_Diversity_PJL': pjl_stats.Haplotype_Diversity
                })
            else:
                snp.update({
                    'TajimaD_PJL': 'N/A',
                    'Nucleotide_Diversity_PJL': 'N/A',
                    'Haplotype_Diversity_PJL': 'N/A'
                })

            # Fetch comparative statistics
            comp_query = text("""
                SELECT
                    FST_Score,
                    XPEHH_Score
                FROM Comparative_Statistics
                WHERE RSID = :snp_id
            """)
            comp_stats = conn.execute(comp_query, {'snp_id': snp_id}).fetchone()
            if comp_stats:
                snp.update({
                    'FST_Score': comp_stats.FST_Score,
                    'XPEHH_Score': comp_stats.XPEHH_Score
                })
            else:
                snp.update({
                    'FST_Score': 'N/A',
                    'XPEHH_Score': 'N/A'
                })

            return render_template("snp_detail.html", snp=snp)
            
    except Exception as e:
        logger.error(f"Error fetching SNP details: {e}")
        return render_template("snp_detail.html", snp=None)
    

@app.route("/download_snp_data", methods=["GET"])
def download_snp_data():
    # gets search items
    search_type = request.args.get("filter")
    query = request.args.get("query")
    chromosome = request.args.get("chromosome")
    start_position = request.args.get("start_position")
    end_position = request.args.get("end_position")
    populations = request.args.getlist("Population")
    phenotypes = request.args.getlist("phenotype")

    # runs execute query to get results
    total_results, results = execute_query(search_type, query, populations, phenotypes, chromosome, start_position, end_position, page=1, per_page=1000000)

    # puts results collected into a df via pandas
    df = pd.DataFrame(results)

    # makes text file with IO import
    output = io.StringIO()
    df.to_csv(output, index=False, sep='\t')
    output.seek(0)

    # Create a response with the text file
    response = Response(
        output,
        mimetype="text/plain",
        headers={"Content-Disposition": f"attachment;filename=snp_data.txt"}
    )

    return response

# About route
@app.route("/about")
def about():
    return render_template("about.html")

# User Guide route
@app.route("/userguide")
def userguide():
    return render_template("userguide.html")

# Contact Us route
@app.route("/contactus")
def contact():
    return render_template("contact.html")

# Data Sources route
@app.route("/datasources")
def datasources():
    return render_template("datasources.html")

@app.route("/survey", methods=["GET", "POST"])
def survey():
    if request.method == "POST":
        result = request.form
        return render_template("survey.html", result=result)

    return render_template("survey.html", result=None)

# KEGG API Functions
def get_kegg_gene_id(gene_name):
    """
    Use KEGG's API to convert a gene symbol to a KEGG gene ID e.g., for gene "TCF7L2", will be "hsa:6934".
    """
    kegg_find_api = f"http://rest.kegg.jp/find/hsa/{gene_name}"  # KEGG'S API endpoint using gene symbol
    try:
        response = requests.get(kegg_find_api, timeout=5).text  # retrieves API response as text
        logger.info(f"KEGG find API response for {gene_name}: {response}")  # log the API response for easier debugging
        if response:  # ensure the response is a valid one
            first_line = response.splitlines()[0]  # extracts KEGG gene ID e.g. "hsa:6934"
            kegg_gene_id = first_line.split("\t")[0].strip()  # e.g., "hsa:6934"
            logger.info(f"Extracted KEGG gene ID for {gene_name}: {kegg_gene_id}")  # log the extracted gene ID for easier debugging
            return kegg_gene_id
    except Exception as e:
        logger.error(f"KEGG find API error: {e}")  # catches and logs API errors for easier debugging
        return None  # returns nothing if no ID is found or an error occurs

# Fetch gene functional data through KEGG's API
@app.route("/fetch_gene_data")
def fetch_gene_data():
    """
    Redirects user to a new page when clicking on a gene - user will be redirected to KEGG's gene entry page for functional information on that gene.
    """
    gene_name = request.args.get("gene")  # get gene name for the request's query parameters
    logger.info(f"Fetching gene data for: {gene_name}")  # logs the gene name that is being processed for easier debugging
    # convert gene symbol to KEGG gene ID using get_kegg_gene_id function
    kegg_gene_id = get_kegg_gene_id(gene_name)
    # only proceeds if a valid KEGG ID is found
    if kegg_gene_id:
        logger.info(f"Using KEGG gene ID: {kegg_gene_id}")  # logs the KEGG ID that has been found for easier debugging
        return redirect(f"https://www.kegg.jp/entry/{kegg_gene_id}")  # redirects user to new page with KEGG gene entry
    # fall back to Ensembl GRCh37 gene summary page if no KEGG data is available
    logger.info(f"No KEGG gene ID found for {gene_name}, redirecting to Ensembl GRCh37 for {gene_name}")
    return redirect(f"https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_name}")

# Run the Flask app
if __name__ == "__main__":
    if test_db_connection():
        app.run(debug=True)
    else:
        logger.error("Failed to connect to the database. Exiting.")