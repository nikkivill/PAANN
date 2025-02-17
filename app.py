from flask import Flask, render_template, request
import json


app = Flask(__name__)

with open("SNP_table.json", "r") as file:
    snp_data = json.load(file)

@app.route("/")
def home():
    return render_template("home.html")

def search_snp(query, search_type, populations):
    matches = []
    
    for snp in snp_data:
        # Filter by search query
        if search_type == "snp_id" and query.lower() in snp["RSID"].lower():
            matches.append(snp)
        elif search_type == "chromosome" and str(query) == str(snp["Chromosome"]):
            matches.append(snp)
        elif search_type == "position" and str(query) == str(snp["Position"]):
            matches.append(snp)

    # Apply population filter
    if populations:
        filtered_matches = [snp for snp in matches if snp["Selection_Population"] in populations]
        return filtered_matches
    return matches


@app.route("/search", methods=["GET"])
def search():
    query = request.args.get("query")
    search_type = request.args.get("filter")
    populations = request.args.getlist("Population")  # Get selected populations

    if query and search_type:
        results = search_snp(query, search_type, populations)
        return render_template("search_results.html", query=query, results=results, search_type=search_type)

    return render_template("search_results.html", query=None, results=[], search_type=None)


@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/userguide")
def userguide():
    return render_template("userguide.html")

@app.route("/contactus")
def contact():
    return render_template("contact.html")

@app.route("/datasources")
def datasources():
    return render_template("datasources.html")

@app.route("/survey", methods=["GET", "POST"])
def survey():
    if request.method == "POST":
        age = request.form.get("age")
        family_history = request.form.get("family_history")
        bmi = request.form.get("bmi")
        physical_activity = request.form.get("physical_activity")
        smoking = request.form.get("smoking")
        diet = request.form.get("diet")

        # Basic Risk Assessment Logic
        risk_factors = 0
        if int(age) > 45:
            risk_factors += 1
        if family_history == "yes":
            risk_factors += 2
        if float(bmi) > 30:
            risk_factors += 2
        if physical_activity == "low":
            risk_factors += 1
        if smoking == "yes":
            risk_factors += 1
        if diet == "unhealthy":
            risk_factors += 1

        # Generate Result
        if risk_factors >= 5:
            result = "High Risk for Heart Disease & Hypertension"
        elif risk_factors >= 3:
            result = "Moderate Risk for Metabolic Syndrome"
        else:
            result = "Low Risk, Maintain a Healthy Lifestyle"

        return render_template("survey.html", result=result)

    return render_template("survey.html", result=None)

if __name__ == "__main__":
    app.run(debug=True)