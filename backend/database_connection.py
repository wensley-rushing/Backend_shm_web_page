import sys
from pathlib import Path

import pandas as pd
import psycopg2
import yaml
from scipy import stats, signal, sparse
from sqlalchemy import create_engine,exc

## Global Settigns and credentials loads
current_path = Path(str(__file__)).parent

with open(current_path/'credentials.yaml') as file:
    credentials = yaml.load(file, Loader=yaml.FullLoader)
    

def preprocess(df,columns):
    """
    Applies detrend to all rows in sensors df. 
    Parameters
    ----------
    df: DataFrame
        Dataframe containing sensors data.
    
    columns: List[str]
        List containing the name of the columns to use.
        
    Returns
    -------
    _df: DataFrame
        Returns Dataframe containing detrended dataset 
    """
    _df = df[columns].copy()
    _df[columns[1:]] = df[columns[1:]].apply(signal.detrend)
    return _df

def create_connection():
    """
    Creates the engine to connect to the Postgres DB
    
    Parameters
    ----------
    None
        
    Returns
    -------
    engine: sqlalchemy object
        Returns the engine to connect to the database.
    """
    engine = None
    
    host = credentials['host']
    database = credentials['database_name']
    user = credentials['user']
    password = credentials['password']
    port = credentials['port']
          
    try:
        engine = create_engine(f'postgresql://{user}:{password}@{host}:{port}/{database}')
    except exc.SQLAlchemyError as e:
        print(f'Error {e}')
        sys.exit(1)
        
    return engine

def connection_and_data_retrieving(Q1,Q2):
    """
    Stablish a connection with the database and returns the data
    between Q1 and Q2.
    
    Parameters
    ----------
    Q1, Q2: datetime 
    
    Returns
    -------
    accelerations_df: DataFrame
        Returns the DataFrame with the data between Q1 and Q2
    create_connection: Function
        Returns the engine that connects to the Postgres DB.
    """
    con = None
    try:
        con = psycopg2.connect(
            host=credentials['host'],
            database=credentials['database_name'],
            user=credentials['user'],
            password=credentials['password']  
        )
    except psycopg2.DatabaseError as e:
        print(f'Error {e}')
        sys.exit(1)

    df = pd.read_sql(
    f'''
    SELECT *
    FROM "100hz_python_SIM" sim
    WHERE sim.sampledatetime
        BETWEEN '{Q1}'
            AND '{Q2}';
    ''', con)
    columns = ['sampledatetime','AI1','AI2','AI3','AI4','AI5','AI10','AI11','AI12']
    accelerations_df = preprocess(df,columns)
    con.close()
    
    return accelerations_df,create_connection()