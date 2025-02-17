from flask import Flask, render_template, request
import pandas as pd
from sqlalchemy import create_engine, text
import logging
import math


app = Flask(__name__)

# Logging allows us to track events that happen when flask runs.
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# MYSQL Database connection : P2
DB_USERNAME = "root"
DB_PASSWORD = "Family786"
DB_HOST = "localhost"
DB_NAME = "P2"

# Create a database engine to connect to the database server : P2, did f string so easier changes in future
try:
    engine = create_engine(
        f"mysql+mysqlconnector://{DB_USERNAME}:{DB_PASSWORD}@{DB_HOST}/{DB_NAME}",
        pool_pre_ping=True
    )
    logger.info("Database engine created successfully.")
except Exception as e:
    logger.error(f"Error creating database engine: {e}")
    engine = None

# Function to test the database connection before running the app
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
    # checks if request is POST, gets query from form data, gets dropdown filter value, and gets what population, same for GET
    if request.method == "POST":
        query = request.form.get("query")
        search_type = request.form.get("filter")
        populations = request.form.getlist("Population")
        phenotypes = request.form.getlist("phenotype")
    else: 
        query = request.args.get("query")
        search_type = request.args.get("filter")
        populations = request.args.getlist("Population")
        phenotypes = request.args.getlist("phenotype")

    # gets a page parameter, the current one, and sets how many results per page
    page = request.args.get('page', 1, type=int)
    per_page = 20

    # if both query and search type given, executes search query, by calling on execute query() function to do the aa database search 
    # Also does pagnation so calculates total pages with maths import, and returns the search results
    if query and search_type:

        search_value = query + "%"
        total_results, results = execute_query(search_type, search_value, populations, phenotypes, page, per_page)
        total_pages = math.ceil(total_results / per_page)
        logger.info(f"Results: {results}")

        return render_template("search_results.html", query=query, results=results, search_type=search_type, page=page, total_pages=total_pages)

    return render_template("search_results.html", query=None, results=[], search_type=None, page=1, total_pages=1)

def execute_query(search_type, search_value, populations=None, phenotypes=None, page=1, per_page=20):
    # Dictionary to map search types to SQL, LIKE used for partial matches so rs... gives return rs1234, rs5678, etc.
    query_map = {
        "snp_id": "sc.RSID LIKE :search_value",  # Use LIKE for partial matches
        "gene": "g.Mapped_Gene LIKE :search_value",
        "position": "sc.Position LIKE :search_value",
        "chromosome": "sc.Chromosome LIKE :search_value"
    }

    if search_type not in query_map:
        logger.error(f"Invalid search type: {search_type}")
        return 0, []

    # Base SQL query!!! joins all those tables, and selects the columns we want, based on the search type and search map
    sql_query = f"""
        SELECT
            sc.SNP_ID,
            sc.RSID,
            sc.Chromosome,
            sc.Position,
            sc.Reference_Allele,
            sc.Alternate_Allele,
            sc.P_value,
            g.Mapped_Gene,
            ph.Phenotype AS Phenotype,
           COALESCE(p.Population_ID, 'N/A') AS Population
        FROM SNP sc
        LEFT JOIN Gene g ON sc.Gene_ID = g.Gene_ID
        LEFT JOIN Phenotype ph ON sc.SNP_ID = ph.SNP_ID
        LEFT JOIN SNP_Population sp ON sc.SNP_ID = sp.SNP_ID
        LEFT JOIN Population p ON sp.Population_ID = p.Population_ID
        WHERE {query_map[search_type]}
    """
    
    # Parameters for the query
    params = {"search_value": f"%{search_value}%"}

    # if populations selected, filters them, converts into tuple for better query
    if populations:
        placeholders = ",".join([f":pop{i}" for i in range(len(populations))])
        sql_query += f" AND sp.Population_ID IN ({placeholders})"
    for i, pop in enumerate(populations):
        params[f"pop{i}"] = pop


    # if phenotypes selected, filters them, converts into tuple for better query
    if phenotypes:
        sql_query += f" AND ph.Phenotype IN ({','.join([':ph' + str(i) for i in range(len(phenotypes))])})"
    for i, ph in enumerate(phenotypes):
        params[f"ph{i}"] = ph

    # This part is ai coded so gotta understand it and change it later, but basically goes the pages
    try:
        with engine.connect() as conn:
            logger.info(f"Executing query: {sql_query}")
            logger.info(f"Query parameters: {params}")

            # First, get the total number of results
            count_query = f"SELECT COUNT(*) FROM ({sql_query}) AS subquery"
            total_results = pd.read_sql(text(count_query), conn, params=params).iloc[0, 0]

            # Modify the query to include pagination
            sql_query += " LIMIT :limit OFFSET :offset"
            params["limit"] = per_page
            params["offset"] = (page - 1) * per_page

            # Execute the paginated query and fetch results
            results = pd.read_sql(text(sql_query), conn, params=params)
            logger.info(f"Query returned {len(results)} rows.")

            return total_results, results.to_dict(orient='records')
    except Exception as e:
        logger.error(f"Error executing query: {e}")
        return 0, []  # Return 0 total results and an empty list

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

# Survey route
@app.route("/survey", methods=["GET", "POST"])
def survey():
    if request.method == "POST":
        age = request.form.get("age")
        family_history = request.form.get("family_history")
        bmi = request.form.get("bmi")
        physical_activity = request.form.get("physical_activity")
        smoking = request.form.get("smoking")
        diet = request.form.get("diet")

        # Calculate risk factors
        risk_factors = sum([
            int(age) > 45,
            family_history == "yes",
            float(bmi) > 30,
            physical_activity == "low",
            smoking == "yes",
            diet == "unhealthy"
        ])

        # Determine risk level
        result = (
            "High Risk for Heart Disease & Hypertension" if risk_factors >= 5 else
            "Moderate Risk for Metabolic Syndrome" if risk_factors >= 3 else
            "Low Risk, Maintain a Healthy Lifestyle"
        )

        return render_template("survey.html", result=result)

    return render_template("survey.html", result=None)

# Run the Flask app
if __name__ == "__main__":
    if test_db_connection():
        app.run(debug=True)
    else:
        logger.error("Failed to connect to the database. Exiting.")
