import pandas as pd

def map_data(snp_file, inserting_file, output_file, columns_to_drop):
    """This function allows you to map your CSV file e.g. statistic files, to the correct SNP_ID
    and allow you to drop unnecessary columns before intrgating data into the schema.
    
    The parameters are:
    - snp_file: this is the exported SNP_mapping.csv file from the schema that includes the SNP_ID fields
    - inserting_file: this is the csv file of data that needs to be inserted 
    - output_file: this is the name of the new file once merged 
    - columns_to_drop: these are the list of column names that can be dropped in the new file
    """
    #loading the two files to merge
    snp = pd.read_csv(snp_file)
    inserting_file = pd.read_csv(inserting_file)

    #mapping the inserting file into snp mapping csv file on 'rsid' field 
    mapped = snp.merge(inserting_file, on='RSID', how='left')

    # dropping extra columns not needed
    mapped = mapped.drop(columns=columns_to_drop)
    
    # saving the merged csv into a new output file
    mapped.to_csv(output_file, index=False)

