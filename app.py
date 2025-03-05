from flask import Flask, render_template, request, Response, redirect, flash, url_for  # Importing Flask and other packages like redirects, flash messages, etc
from flask_sqlalchemy import SQLAlchemy # Importing SQLAlchemy for database management
from sqlalchemy import create_engine, text # Importing create_engine and text from SQLAlchemy
import pandas as pd # Importing pandas for data manipulation
import logging # Importing logging for debugging
import math # Importing math for calculations
import io # Importing io for file handling
import requests # Importing requests for HTTP requests
from datetime import datetime  # Importing datetime for timestamps

# Sets up logging for debugging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize Flask app
app = Flask(__name__)

# Set a secret key for session management
app.secret_key = 'secret_key' 

# a function to create a database engine
def create_db_engine(user, password, host, database):
    return create_engine(f"mysql+mysqlconnector://{user}:{password}@{host}/{database}")

# Prompts the user for database details
user = input("Enter MySQL username: ")
password = input("Enter MySQL password: ")
host = input("Enter MySQL host: ")
database = input("Enter MySQL database: ")

# Creates the database engine with the above details using previous function
engine = create_db_engine(user, password, host, database)

# Tests database connection and presence of initial engine
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

# Makes sure SQLAlchemy uses MySQL database
app.config['SQLALCHEMY_DATABASE_URI'] = f"mysql+mysqlconnector://{user}:{password}@{host}/{database}"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # This extra step improves performance by removing modification tracking
db = SQLAlchemy(app)

# Defines Message Model for storing user messages in the database
class Message(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False)
    email = db.Column(db.String(100), nullable=False)
    message = db.Column(db.Text, nullable=False)
    created_at = db.Column(db.DateTime, default=datetime.utcnow)

# Create all database tables if they dont already exist
with app.app_context():
    db.create_all()

# Base SQL query template, structured so mul;tiple cna be added later on and stord efficiently for later use
QUERIES = {
    "base_snp_query": """
        SELECT
            s.SNP_ID,
            s.RSID,
            s.Chromosome,
            s.Position,
            s.Reference_Allele,
            s.Alternate_Allele,
            CONCAT(ROUND(s.P_value / POW(10, FLOOR(LOG10(s.P_value))), 2), 'e', FLOOR(LOG10(s.P_value))) AS P_value,
            s.Mapped_genes,
            s.Phenotype,
            COALESCE(st.Population_ID, 'N/A') AS Population,
            cs.FST_Score,
            cs.XPEHH,
            st.Tajima_D,
            st.Haplotype_Diversity,
            st.Nucleotide_Diversity,
            s.Ancestry
        FROM SNP_Information s
        LEFT JOIN Statistics st ON s.SNP_ID = st.SNP_ID
        LEFT JOIN Population p ON st.Population_ID = p.Population_ID
        LEFT JOIN Comparative_Statistics cs ON s.SNP_ID = cs.SNP_ID
    """
}

# Build query function that builds query based on search parameters (e.g., genomic region, gene name)
def build_query(search_type, search_value, populations=None, phenotypes=None, 
                chromosome=None, start_position=None, end_position=None):
    """Build SQL query based on search parameters."""
    query = QUERIES["base_snp_query"] + " WHERE " # Gets base query and appends "WHERE" to start adding filtering conditions based on user's search parameters
    params = {} # Creates an empty dictionary to store query parameters

    # If the search is a genomic region, query is formatted to inclyde chromosome and positions from MySQL
    if search_type == "genomic_region":
        query += "s.Chromosome = :chromosome AND s.Position BETWEEN :start_position AND :end_position" # Appends SQL condiiton to filter SNP based on chorosome and given positions
        params.update({   # Updates the params dictionary with values needed for the SQL query
            "chromosome": chromosome,
            "start_position": int(start_position) if start_position else 0,
            "end_position": int(end_position) if end_position else 0
        })
    else: # If not genomic region, so else, then return to normal search conditions
        search_conditions = {
            "snp_id": "s.RSID LIKE :search_value",
            "gene": "s.Mapped_genes LIKE :search_value",
            "position": "s.Position LIKE :search_value",
            "chromosome": "s.Chromosome LIKE :search_value"
        }
        query += search_conditions.get(search_type, "1=1") # Appends SQL condition based on the search type
        params["search_value"] = f"%{search_value}%" # Adds search value to the params dictionary

    # Add population and phenotype filters by changing the query with placeholders
    if populations:
        pop_placeholders = ",".join([f":pop{i}" for i in range(len(populations))]) # List comprehension to create list of placeholders
        query += f" AND st.Population_ID IN ({pop_placeholders})" # Adds an SQL condition - filters results based on Pop ID
        params.update({f"pop{i}": pop for i, pop in enumerate(populations)}) # Maps each pop to placeholder so : pop0, :BEB, etc

    if phenotypes: # repeat for phenotypes
        ph_placeholders = ",".join([f":ph{i}" for i in range(len(phenotypes))])
        query += f" AND s.Phenotype IN ({ph_placeholders})"
        params.update({f"ph{i}": ph for i, ph in enumerate(phenotypes)})

    return query, params # Returns the final query and params dictionary

# Function to execute query, this takes into account paged results as it is not fetching all results
def execute_query(query, params, fetch_all=False, page=1, per_page=20):

    if not engine:
        logger.error("Database engine not available")
        return 0, []

    try:
        with engine.connect() as conn:
            # Counts the total results for pagination via wrapping the query 
            count_query = f"SELECT COUNT(*) FROM ({query}) AS subquery"
            total_results = pd.read_sql(text(count_query), conn, params=params).iloc[0, 0] # pandas to execute count query, extracts first row, first column so the total count

            # If fetch all is set to false, which it would be, then it limits how many results based on results per page
            if not fetch_all:
                query += " LIMIT :limit OFFSET :offset" # calculates starting row for the page
                params.update({"limit": per_page, "offset": (page - 1) * per_page})

            # Executes query using pandas to read and return via a pandas dataframe
            results = pd.read_sql(text(query), conn, params=params)
            return total_results, results.to_dict(orient='records')
        
    except Exception as e:
        logger.error(f"Error executing query: {e}")
        return 0, []

# Function that fetches and combines sum statistics and comaprative stats
def process_result_with_stats(connection, result):

    # SQL queries to fetch sum stats and comparative stats
    stats_query = """
        SELECT Tajima_D, Nucleotide_Diversity, Haplotype_Diversity 
        FROM Statistics WHERE SNP_ID = :snp_id AND Population_ID = :population
    """
    comparative_query = """
        SELECT FST_Score, XPEHH 
        FROM Comparative_Statistics WHERE SNP_ID = :snp_id
    """

    # Fetches stats for PJL and BEB by executing stats_query fpr PJL and BEB, uses fetchone to get the first row
    summary_pjl = connection.execute(text(stats_query), {'snp_id': result['SNP_ID'], 'population': 'PJL'}).fetchone()
    summary_beb = connection.execute(text(stats_query), {'snp_id': result['SNP_ID'], 'population': 'BEB'}).fetchone()
    comparative_stats = connection.execute(text(comparative_query), {'snp_id': result['SNP_ID']}).fetchone()

    # Combines results into a new dictionary to be called on later, uses if else statemnts if none values returned
    combined_result = {
        **result,
        'TajimaD_PJL': summary_pjl.Tajima_D if summary_pjl else None,
        'Nucleotide_Diversity_PJL': summary_pjl.Nucleotide_Diversity if summary_pjl else None,
        'Haplotype_Diversity_PJL': summary_pjl.Haplotype_Diversity if summary_pjl else None,
        'TajimaD_BEB': summary_beb.Tajima_D if summary_beb else None,
        'Nucleotide_Diversity_BEB': summary_beb.Nucleotide_Diversity if summary_beb else None,
        'Haplotype_Diversity_BEB': summary_beb.Haplotype_Diversity if summary_beb else None,
        'FST_Score': comparative_stats.FST_Score if comparative_stats else None,
        'XPEHH': comparative_stats.XPEHH if comparative_stats else None,
    }
    return combined_result

# Function to get the search parameters of the user
def get_search_params():

    form = request.form if request.method == "POST" else request.args
    search_params = {
        "search_type": form.get("filter"), # field to filter for
        "search_value": form.get("query"), # value to search for
        "populations": form.getlist("Population") or ["BEB", "PJL"],
        "phenotypes": form.getlist("phenotype") or ["FGadjBMI", "HBA1CadjBMI", "FIadjBMI", "2hrGadjBMI", "T2DadjBMI"],
        "chromosome": form.get("chromosome"),
        "start_position": form.get("start_position"),
        "end_position": form.get("end_position"),
        "page": request.args.get("page", 1, type=int)
    }
    logger.info(f"Search parameters: {search_params}")  # Log the search parameters
    return search_params

# Function to get uniprot accession id via gene name
def get_uniprot_accession_id(gene_name):
# constrcturs UniProt URL with gene name
    uniprot_api = f"https://rest.uniprot.org/uniprotkb/search?query={gene_name}&fields=accession&format=json"

    try:
        response = requests.get(uniprot_api, timeout=5)  # sends request to the uniprot API
        response.raise_for_status()  # if request is not succesful, raise an error and do not proceed
        data = response.json()  # converts response from json to python dictionary

        # makes sure theres a result and extract the first uniprot accession ID
        if "results" in data and data["results"]:
            uniprot_id = data["results"][0]["primaryAccession"]
            return uniprot_id
    # catch and log errors 
    except requests.exceptions.RequestException as e:
        logger.error(f"UniProt API error for {gene_name}: {e}")  

    return None  # return nothing if no valid ID is found or an error occurs

# Initial landing page for site using home.html
@app.route("/")
def home():
    return render_template("home.html")

# Search route with GET and POST, uses the search parameters given to check if valid results returnable
@app.route("/search", methods=["GET", "POST"])
def search():
    # gets and validates parameters
    search_params = get_search_params()
    if search_params["search_type"] != "genomic_region" and not search_params["search_value"]: # if btoh empty, handles with an error message
        return render_template(
            "search_results.html",
            error=f"Please provide a valid search term for {search_params['search_type']} search.",
            **search_params,
            results=[],
            all_results=[],
            total_pages=0
        )  # returns error f' string if query is not in MYSQL database

    # Builds and executes queries based on users search parameters
    query, params = build_query(
        search_type=search_params["search_type"],
        search_value=search_params["search_value"],
        populations=search_params["populations"],
        phenotypes=search_params["phenotypes"],
        chromosome=search_params["chromosome"],
        start_position=search_params["start_position"],
        end_position=search_params["end_position"]
    )

    # Executes query function, but also takes into account pageinated results which is used later for the table view
    total_results, paginated_results = execute_query(query, params, page=search_params["page"], per_page=20)

    # Executes query function, but fecthes all so thta all results are fetched for later use when displaying, for example, all sum stats from user query
    _, all_results = execute_query(query, params, fetch_all=True)

    # Processes results with stats
    with engine.connect() as conn:
        paginated_results = [process_result_with_stats(conn, result) for result in paginated_results]
        all_results = [process_result_with_stats(conn, result) for result in all_results]
    
    # ensures correct rounding if decimal value so all results shown
    total_pages = math.ceil(total_results /20) if total_results > 0 else 0
    return render_template(
        "search_results.html",
        query=search_params["search_value"],
        results=paginated_results,
        all_results=all_results,
        total_results=total_results,
        total_pages=total_pages,
        **search_params
    )

# Function to get snp details via snp ID and population BEB/PJL
def get_snp_details(snp_id, population):
    """Gets detailed information for a specific SNP."""
    if not engine:
        logger.error("Database engine not available")
        return None

    try:
        with engine.connect() as conn:
            # Uses the base query but also checks for population ID and RSID with  a WHERE clause 
            query = QUERIES["base_snp_query"] + " WHERE s.RSID = :snp_id AND st.Population_ID = :population"
            logger.info(f"Executing query: {query} with snp_id={snp_id}, population={population}")  # Debugging log
            result = conn.execute(text(query), {"snp_id": snp_id, "population": population}).fetchone()
            if result:
                return dict(result._mapping)
            else:
                logger.warning(f"No data found for SNP ID: {snp_id}, Population: {population}")
                return None
            
    except Exception as e:
        logger.error(f"Error fetching SNP details: {e}")
        return None
    
# App route for snp details page when clicking 'view' button
@app.route("/snp/<snp_id>")
def snp_detail(snp_id):
    """Route to display detailed information about a specific SNP."""
    population = request.args.get("population") # This page is based on population row (BEB/PJL)
    if not population:
        logger.error("Population parameter is missing.")
        return render_template("snp_detail.html", snp=None)

    # Fetches SNP details, handles errors
    snp = get_snp_details(snp_id, population)
    if not snp:
        logger.warning(f"No SNP found for RSID: {snp_id}, Population: {population}")
        return render_template("snp_detail.html", snp=None)

    return render_template("snp_detail.html", snp=snp)
    
# App route for downloading SNP data using get from the parameters users used
@app.route("/download_snp_data", methods=["GET"])
def download_snp_data():
    # Gets search parameters directly from the users request
    search_type = request.args.get("filter")
    search_value = request.args.get("query")
    chromosome = request.args.get("chromosome")
    start_position = request.args.get("start_position")
    end_position = request.args.get("end_position")
    populations = request.args.getlist("Population")
    phenotypes = request.args.getlist("phenotype")
    
    # Builds query with the previous build_query function, confirms them
    query, params = build_query(
        search_type=search_type,
        search_value=search_value,
        populations=populations,
        phenotypes=phenotypes,
        chromosome=chromosome,
        start_position=start_position,
        end_position=end_position
    )
    
    # Executse query again to get all results (fetch_all=True)
    total_results, results = execute_query(query, params, fetch_all=True)
    
    # Processes results with stats
    with engine.connect() as conn:
        results = [process_result_with_stats(conn, result) for result in results]
    
    # Converts these results to DataFrame and uses StringIO, a package used to handle input/output files
    df = pd.DataFrame(results)
    output = io.StringIO()
    
    # IF else statement to handle empty pd.df
    if not df.empty:
        # Finds unique SNPs for statistics calculation
        unique_snps = df.drop_duplicates(subset=['RSID'])
        
        # Calculates average statistics and standard deviations working with unique SNPs, so duplicates nto taken into account
        avg_tajima_d = unique_snps['Tajima_D'].mean() if 'Tajima_D' in unique_snps.columns else None
        sd_tajima_d = unique_snps['Tajima_D'].std() if 'Tajima_D' in unique_snps.columns else None
        
        avg_nuc_div = unique_snps['Nucleotide_Diversity'].mean() if 'Nucleotide_Diversity' in unique_snps.columns else None
        sd_nuc_div = unique_snps['Nucleotide_Diversity'].std() if 'Nucleotide_Diversity' in unique_snps.columns else None
        
        avg_hap_div = unique_snps['Haplotype_Diversity'].mean() if 'Haplotype_Diversity' in unique_snps.columns else None
        sd_hap_div = unique_snps['Haplotype_Diversity'].std() if 'Haplotype_Diversity' in unique_snps.columns else None
        
        avg_fst = unique_snps['FST_Score'].mean() if 'FST_Score' in unique_snps.columns else None
        sd_fst = unique_snps['FST_Score'].std() if 'FST_Score' in unique_snps.columns else None
        
        avg_xpehh = unique_snps['XPEHH'].mean() if 'XPEHH' in unique_snps.columns else None
        sd_xpehh = unique_snps['XPEHH'].std() if 'XPEHH' in unique_snps.columns else None
        
        # Creates the  headers for region statistics
        stats_headers = [
            "Region_Type", "Chromosome", "Start_Position", "End_Position", "SNP_Count", 
            "Avg_TajimaD", "SD_TajimaD", 
            "Avg_Nucleotide_Diversity", "SD_Nucleotide_Diversity", 
            "Avg_Haplotype_Diversity", "SD_Haplotype_Diversity", 
            "Avg_FST_Score", "SD_FST_Score", 
            "Avg_XPEHH_Score", "SD_XPEHH_Score"
        ]
        
        # Registers the values to be lined up with the headers
        stats_values = [
            "Region_Statistics", chromosome, start_position, end_position, len(unique_snps),
            avg_tajima_d, sd_tajima_d,
            avg_nuc_div, sd_nuc_div,
            avg_hap_div, sd_hap_div,
            avg_fst, sd_fst,
            avg_xpehh, sd_xpehh
        ]
        
        # Writes headers 
        output.write('\t'.join(stats_headers) + '\n')
        
        # Converts None values to "N/A"
        stats_values = ["N/A" if v is None else v for v in stats_values]
        output.write('\t'.join(map(str, stats_values)) + '\n')
        
        # Adds a separator line
        output.write('\n')
        
        # Reorders the columns to move RSID to the first position
        if 'RSID' in df.columns:
            columns_order = ['RSID'] + [col for col in df.columns if col != 'RSID' and col != 'SNP_ID']
            df = df[columns_order]
        
        # Writes SNP data
        df.to_csv(output, index=False, sep='\t', mode='a')
    else:
        # Whilst handling empty results 
        output.write("No results found for the specified criteria\n\n")

    # Ensures pointer starts from top    
    output.seek(0)
    
    # Names the downloaded file based on the searched genomic region
    if chromosome and start_position and end_position:
        filename = f"snp_data_chr{chromosome}_{start_position}_{end_position}.tsv"
    else:
        filename = "snp_data.tsv"
    
    # Returns the file and makes sure its a text file and is named with the previous region file name 
    return Response(
        output.getvalue(),
        mimetype="text/tab-separated-values",
        headers={"Content-Disposition": f"attachment;filename={filename}"}
    )

# Route to fetch gene data
@app.route("/fetch_gene_data")
def fetch_gene_data():
    """
    - Redirects user to a new page when clicking on a gene. 
    - The user is first redirected to PAN-GO entry page for functional information on the clicked gene (GO terms).
    - If no PAN-GO entry is found, it falls back to Ensembl GRCh37.
    """
    gene_name = request.args.get("gene")  # gets gene name from query

    # Retrieves uniprot accession ID for gene name using provided function
    uniprot_id = get_uniprot_accession_id(gene_name)

    # Redirects to PAN-GO if a valid uniprot accession ID is found
    if uniprot_id:
        pan_go_url=f"https://functionome.geneontology.org/gene/UniProtKB:{uniprot_id}"
        response = requests.get(pan_go_url, timeout=5)

        if response.status_code == 200: # ONLY if a PAN-GO entry is found 
            return redirect(pan_go_url) # then redirects to PAN-GO
 
    # otherwise falls back to Ensembl GRCh37 if no PAN-GO entry is found 
    return redirect(f"https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_name}")

# Static pages
@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/userguide")
def userguide():
    return render_template("userguide.html")

# Route to display the contact form
@app.route("/contact", methods=["GET", "POST"])
def contact():
    # User enters so makes this a POST
    if request.method == "POST":
        name = request.form.get("name")
        email = request.form.get("email")
        message = request.form.get("message")

        if name and email and message: # If all given, forms a session and returns a flash message for the user
            new_message = Message(name=name, email=email, message=message)
            db.session.add(new_message)
            db.session.commit()
            flash("Message sent successfully!", "success")
            return redirect(url_for("contact"))

    return render_template("contact.html")

# Route to display messages in the admin panel
@app.route("/admin/messages")
def admin_messages():
    messages = Message.query.order_by(Message.created_at.desc()).all()
    return render_template("admin_messages.html", messages=messages)

# Route to delete a message manually
@app.route("/admin/messages/delete/<int:message_id>", methods=["POST"])
def delete_message(message_id):
    message = Message.query.get(message_id)
    if message:
        db.session.delete(message)
        db.session.commit()
        flash("Message deleted successfully!", "success")
    return redirect(url_for("admin_messages"))


# Static pages
@app.route("/datasources")
def datasources():
    return render_template("datasources.html")

@app.route("/survey", methods=["GET", "POST"])
def survey():
    result = request.form if request.method == "POST" else None
    return render_template("survey.html", result=result)

if __name__ == "__main__":
    if engine:
        app.run(debug=True)
    else:
        logger.error("Failed to initialize database engine. Exiting.")