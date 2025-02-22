import pandas as pd
from sqlalchemy import create_engine

def integration(csv_file, table_name, columns=None, engine=None):
    """
    This function takes the data from the csv file and inserts it into the SQL table based on the corresponding field names/columns
    - csv_file: the name of the csv that has the columns to be inserted
    - table_name: the name of the table in MySQL that you want to insert into
    - columns: a list of the fields from the csv file that are being inserted into SQL 
    - engine = database connection - (mysql+mysqlconnector://username:password@host/database name)
    """
    try:
        data = pd.read_csv(csv_file)
        if columns:
            data = data[columns]

        data.to_sql(table_name, engine, if_exists='append', index=False)
        print(f"Data successfully inserted into MySQL {table_name}")
    except Exception as e:
        print(f"Error inserting data into MySQL {table_name}: {e}")
