import pandas as pd
from sqlalchemy import create_engine

def create_db_engine(user, password, host, database):
    """
    Creates and returns a SQLAlchemy engine for MySQL.
    """
    return create_engine(f"mysql+mysqlconnector://{user}:{password}@{host}/{database}")

def integrate_schema(csv_file, table_name, columns=None, engine=None):
    """
    Insert data from a CSV file into a specific table in MySQL.

    Parameters:
    - csv_file: the name of the CSV file to get the data from.
    - table_name: the name of the SQL table where data will be inserted.
    - columns: list of fields to be selected from the CSV before inserting.
    - engine: SQLAlchemy database engine for connection.
    """
    try:
        data = pd.read_csv(csv_file)
        data = data.fillna('NULL')
        
        if columns:
            data = data[columns]

        if engine is None:
            raise ValueError("Database engine not provided.")

        data.to_sql(table_name, engine, if_exists='append', index=False)
        print(f"Data successfully transferred to SQL table {table_name}")

    except Exception as e:
        print(f"Error inserting data into {table_name}: {e}")

user = input("Enter MySQL username: ")
password = input("Enter MySQL password: ")
host = input("Enter MySQL host: ")
database = input("Enter MySQL database: ")

engine = create_db_engine(user, password, host, database)

