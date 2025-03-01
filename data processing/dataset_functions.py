
 ### Function for processing and combining association data CSV files from T2D Knowledge Portal and generating TXT file for rsIDs

import pandas as pd
import os

def process_association_data(input_folder: str, output_csv_name="processed_association_data.csv", output_rsids_name="all_rsids.txt"):
    """
    Processes a folder of multiple CSV files of association data downloaded from T2D Knowledge Portal and combines
    them into a single CSV file.
    
    - Removes rows with missing rsIDs in the 'dbSNP' column
    - Removes rows where chromosome is X or Y.
    - Keeps only specific columns.
    - Orders the data by chromosome and position.
    - Combines processed data into one CSV file (used as the main SNP table)
    - Generates a separate text file of all rsIDs with duplicates removed (for filtering for population specific SNP
    data for statistical tests). 

    Parameters: 
    input_folder (str): Path to the folder containing the CSV files.
    output_csv_name (str, optional): The name of the output CSV file (default: "processed_association_data.csv").
    output_rsids_name (str, optional): The name of the output text file containing the rsIDs (default: "all_rsids.txt").
    
    """
    # empty list to store data from all files
    combined_data=[]
    
    # iterate over all the CSV files within the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".csv"):
            file_path = os.path.join(input_folder, file_name) # find file path 
            
            # read the current csv file
            df = pd.read_csv(file_path)
            
            # remove rows where dbSNP (rsID) is missing 
            df = df.dropna(subset=['dbSNP'])
            # remove rows where the chromosome is X or Y 
            df = df[~df['chromosome'].astype(str).str.upper().isin(['X', 'Y'])]
            
            # keep only specific columns
            keep_columns = ['dbSNP', 'chromosome', 'position', 'alt', 'reference', 'phenotype', 'pValue', 'ancestry']
            df = df[keep_columns]
            
            # rename some columns
            df = df.rename(columns={
                'dbSNP' : 'RSID',
                'chromosome': 'Chromosome',
                'position': 'Position',
                'alt' : 'Alternate_Allele',
                'reference' : 'Reference_Allele',
                'phenotype' : 'Phenotype',
                'pValue' : 'P_value',
                'ancestry' : 'Ancestry'
            })
            
            # append the data to the combined_data list
            combined_data.append(df)
            
    # combine combined_data into a single dataframe
    combined_df = pd.concat(combined_data, ignore_index=True)
    
    # ensure chromosome column is integer before sorting
    combined_df['Chromosome'] = combined_df['Chromosome'].astype(int)

    # sort the df by chromosome and position
    combined_df = combined_df.sort_values(by=['Chromosome', 'Position'])
    
    # output file name based on name provided or default
    output_csv = output_csv_name
    output_rsids_txt = output_rsids_name

    # save the combined data to a CSV file
    combined_df.to_csv(output_csv, index=False)

    # print number of SNPs in output CSV file
    num_snps = len(combined_df)
    print(f"Total number of SNPs in {output_csv}: {num_snps}")

    # extract rsids, drop duplicates, and save as a text file
    unique_rsids = combined_df['RSID'].drop_duplicates()
    unique_rsids.to_csv(output_rsids_txt, index=False, header=False)

    # print number of SNPs in the rsIDs text file
    with open(output_rsids_txt, "r") as file:
        rsid_count = sum(1 for line in file if line.strip())
    print(f"Total number of SNPs in {output_rsids_txt}: {rsid_count}")





### Function for removing missing SNPs (after filtering for PJL and BEB SNP data) from combined association data CSV file

import pandas as pd 

def remove_missing_snps(input_csv: str, input_txt: str, output_csv_name="final_association_data.csv", output_txt_name="final_unique_rsids.txt"):
    """
    Removes SNPs from a CSV file where the rsID matches any rsID in a provided text file.
    
    - Reads a CSV file containing an rsID column and a text file with a list of rsIDs to remove (one rsID per line).
    - Filters out the rows from the CSV where the `rsID` is found in the text file.
    - Saves the data into a new CSV file.
    - Saves the remaining rsIDs into a text file.
    
    Parameters:
    input_csv (str): Path to the input CSV file containing SNP association data.
    input_txt (str): Path to the text file containing rsIDs to be removed, one per line.
    output_csv_name (str, optional): Path for saving the cleaned CSV file (default: "final_association_data.csv").
    output_txt_name (str, optional): Path for saving the remaining rsIDs text file (default: "final_unique_rsids.txt").
    
    """
    # load csv file into a dataframe
    df = pd.read_csv(input_csv)
    
    # load the list of rsIDs from the text file to remove (it should be one rsID per line already)
    with open(input_txt, 'r') as file:
        rsids_to_remove = [line.strip() for line in file.readlines()]
        
    # remove rows where are the rsID is in the rsids_to_remove list
    df_cleaned = df[~df['RSID'].isin(rsids_to_remove)]
    
    # number of SNPs in CSV file after removal
    num_snps_csv = len(df_cleaned)
    
    # save to output_csv_name
    df_cleaned.to_csv(output_csv_name, index=False)

    # extract remaining rsIDs to a text file and drop duplicates (needed for mapped genes function)
    remaining_rsids = df_cleaned['RSID'].drop_duplicates().tolist()
    
    with open(output_txt_name, 'w') as output_txt:
        for rsid in remaining_rsids:
            output_txt.write(f"{rsid}\n")
            
    # number of SNPs in txt file after removal
    num_snps_txt = len(remaining_rsids)
    
    # print number of SNPs left in ouput_csv_name and output_txt_name
    print(f"Total SNPs left in {output_csv_name} after removal: {num_snps_csv}")
    print(f"Total SNPs left in {output_txt_name} after removal: {num_snps_txt}")





### Function for fetching mapped gene(s) (if any) for each SNP rsID from Ensembl GRCh37 VEP API

import pandas as pd
import requests
import time

def fetch_mapped_genes(rsid_file, output_file="vep_mapped_genes.csv", batch_size=50, sleep_time=2):
    """
    Fetch mapped gene names for a list of rsIDs using the Ensembl VEP API (GRCh37).
    For each rsID, the function produces one row per gene it maps to (no duplicates), omitting rsIDs with no mapped genes. 

    Parameters:
    - rsid_file (str): Path to the text file containing the rsIDs, one per line.
    - output_file (str, optional): Path for saving the CSV file with mapped genes results (default: "vep_mapped_genes.csv").
    - batch_size (int, optional): Number of rsIDs to process per API request (default: 50).
    - sleep_time (int, optional): Seconds to pause between API requests to avoid rate limits (default: 2).

    Returns:
    - A CSV file with a column for rsIDs and corresponding mapped gene(s).
    """
    # load the rsIDs from the file
    try:
        with open(rsid_file, "r") as f:
            rsid_list = f.read().splitlines()  # read all rsIDs into a list
    except FileNotFoundError: # if the file isn't found, display an error
        print(f"Error: File not found at {rsid_file}. Check the file path.")
        return None

    def get_mapped_genes_vep(rsids):
        """
        Fetch mapped gene names from the Ensembl VEP API for a batch of rsIDs.
        """
        # define the URL for the Ensembl VEP API (GRCh37)
        url = "http://grch37.rest.ensembl.org/vep/human/id"
        headers = {"Content-Type": "application/json"} # header indicating the format

        results = []  # store unique (rsID, gene_name) 

        seen_pairs = set() # store pairs rsID, gene_name pairs 

        try:
            # make a post request to the VEP API with the provided rsIDs
            response = requests.post(url, headers=headers, json={"ids": rsids}, timeout=250)
            if response.status_code == 200: # if the request is successful - status 200
                data = response.json() # parse through JSON response data
                for entry in data: # iterate through each entry in the response
                    rsid = entry.get("id", "unknown") # get rsIDs, or "unknown" if it doesn't exist
                    if rsid != "unknown": # for known rsids only
                        if "transcript_consequences" in entry: # check for transcript consequences (SNP effect on gene transcript)
                          # if gene_symbol is found, return its value - if not, return empty string
                          genes = [conseq.get("gene_symbol", "") for conseq in entry["transcript_consequences"]
                                  if conseq.get("gene_symbol", "").strip()] # ignore empty gene symbols
                          for gene in genes:
                            pair = (rsid, gene)
                            if pair not in seen_pairs: # only add to results if it isn't already there
                              results.append(pair)  # add unique pair to results
                              seen_pairs.add(pair) 
                          
                    
            else: # if API response is not 200, display a warning 
                print(f"Warning: Non-200 status code {response.status_code} for batch: {rsids}")
        except requests.exceptions.RequestException as e: # if there is request exceptions like network errors 
            print(f"Request failed: {e}")

        return results
    
    # process the rsIDs in batches to avoid API overload
    all_results = [] 
    for i in range(0, len(rsid_list), batch_size): # iterate through the rsID list in batches
        batch = rsid_list[i : i + batch_size]  # batch of rsIDs
        batch_results = get_mapped_genes_vep(batch) # get mapped genes for this batch
        all_results.extend(batch_results) # add batch results to list
        time.sleep(sleep_time) # pause for the specificed amount of time to avoid API rate limits

    # convert the list of tuples to a dataframe
    mapped_genes_df = pd.DataFrame(all_results, columns=["rsID", "gene_name"])

    # save the dataframe to a CSV file
    mapped_genes_df.to_csv(output_file, index=False)
    print(f"rsIDs with mapped genes saved as: {output_file}")

    return mapped_genes_df





### Function for adding mapped genes column to combined association data CSV 

import pandas as pd

def add_mapped_genes_to_csv(input_csv, mapped_genes_csv, output_file="final_data_with_genes.csv"):
    """
    Adds a 'Mapped_genes' column to the final association data CSV file, with genes (if more than one) for each rsID
    placed in a comma-separated list, positioned right after the 'ref' column.

    Parameters:
    - input_csv (str): Path to the final association data CSV file.
    - mapped_genes_csv (str): Path to the VEP mapped genes CSV file. 
    - output_file (str, optional): Path to save the final output file (default: "final_data_with_genes.csv")
    """
    
    # load both files
    df_association = pd.read_csv(input_csv)
    df_mapped_genes = pd.read_csv(mapped_genes_csv)

    # group the genes by rsID and join them into a comma-separated string using lambda - x is a group of gene names
    genes_grouped = df_mapped_genes.groupby('rsID')['gene_name'].apply(lambda x: ', '.join(x)).reset_index()

    # rename the column to "RSID"
    genes_grouped = genes_grouped.rename(columns={'rsID' : 'RSID'})
    # rename the column to "Mapped_genes"
    genes_grouped = genes_grouped.rename(columns={'gene_name': 'Mapped_genes'})

    # merge the grouped genes back into the association data by rsID
    merged_df = pd.merge(df_association, genes_grouped, on='RSID', how='left') # keeps all rsIDs in the df_association even if no genes

    # replace NaN values with empty string 
    merged_df['Mapped_genes'] = merged_df['Mapped_genes'].fillna('')

    # reorder the columns to place 'mapped_genes' right after the 'ref' column 
    columns = list(merged_df.columns)
    ref_index = columns.index('Reference_Allele')  # find the index of the ref column
    columns.insert(ref_index + 1, columns.pop(columns.index('Mapped_genes')))  # move mapped genes one position after ref column 

    # reorder the dataframe columns
    merged_df = merged_df[columns]

    # save to CSV file of specified name
    merged_df.to_csv(output_file, index=False)

    print(f"File saved to: {output_file}")




