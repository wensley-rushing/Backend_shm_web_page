from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from sqlalchemy import create_engine

import backend.database_connection as db_conn

engine = db_conn.create_connection()

australian_time = datetime.now() + timedelta(hours = 10)
likelihood_name = f'{australian_time.year}-{australian_time.month:02}-{(australian_time.day-1):02}-likelihood.txt'

code_paths = Path('/home/ec2-user/code_comp_files')
risk_mgmt_paths = Path('/home/ec2-user/risk_management_results')

likelihood = np.loadtxt(code_paths/likelihood_name)

severity_values = np.array([2,2,2,2,3,3,3,0])

def calculate_risk():
    last_data_query = "SELECT pro,monitored_day FROM crack_monitoring ORDER BY monitored_day DESC LIMIT 1"
    last_monitored_crack = pd.read_sql(last_data_query, engine).to_numpy()
    
    pro = last_monitored_crack[0][0]
    date = last_monitored_crack[0][1]
    
    data = np.column_stack((likelihood, pro,severity_values))
    data_df = pd.DataFrame(data, columns = ['ultimate failure','failure due to fatigue', 'severity_values'])
    
    data_df['ultimate failure qualitative'] = pd.cut(
        data_df['ultimate failure'], 
        bins = [-0.1,0.25,0.5,0.75,1.1],
        labels = [
            'Unlikely to happen',
            'Possibly could happen',
            'Likely to happen',
            'Very likely to happen']
    )
    data_df['failure due to fatigue qualitative'] = pd.cut(
        data_df['failure due to fatigue'],
        bins = [-0.1,0.25,0.5,0.75,1.1],
        labels=[
            'Unlikely to happen',
            'Possibly could happen',
            'Likely to happen',
            'Very likely to happen']
    )
    file_text = 'table_data_'+date.strftime('%Y-%m-%d')+'.csv'
    
    data_df.to_csv(risk_mgmt_paths/file_text,index=False)