import math
import os
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from sqlalchemy import create_engine

m = 3.202; #material const of the Paris law
number_of_runs = 2; #number of runs for Monte Carlo Simulation
a_c = [0.35,0.00318,0.00318,0.00318,0.00318,0.00318,0.00318,0.00318];     #rack critical size
Y = 1.12;       #Y constant
C_m = 1.03e-10; #mean of C material variable
C_COV = 0.2;    #covariance of C material variable
C_mu = math.log(C_m**2/(math.sqrt((C_m*C_COV)**2+C_m**2)));  #mean of logarithmic values C material variable
C_sigma = math.sqrt(math.log((((C_m*C_COV)**2)/C_m**2)+1));  #standard deviation of logarithmic values C material variable

date_yesterday_start = (datetime.now() - timedelta(days = 1)).strftime('%Y-%m-%d 00:00:00')
date_yesterday_end = (datetime.now() - timedelta(days = 1)).strftime('%Y-%m-%d 23:59:59')

P_cracks = list()
a_mu_cracks = list()
a_sigma_cracks = list()

def load_data(engine):
    rainflow_data = pd.read_sql(f'''
                            SELECT * 
                            FROM stress_events_rainflow 
                            WHERE timestamp
                            BETWEEN '{date_yesterday_start}'
                            AND '{date_yesterday_end}';
                            ''', engine)
    rainflow_data.drop(columns = ['timestamp'], inplace = True)
    cycles_sum = rainflow_data.sum(axis = 0)
    equivalent_range_mean = rainflow_data.mean(axis = 0)
    equivalent_range_std = rainflow_data.std(axis = 0)
    cycles_sum = cycles_sum.to_numpy()[0::2]
    equivalent_range_mean = equivalent_range_mean.to_numpy()[1::2]
    equivalent_range_std = equivalent_range_std.to_numpy()[1::2]
    
    crack_yesterday = pd.read_sql(f'''
                        SELECT mean_crack, std_crack_size 
                        FROM crack_monitoring 
                        WHERE monitored_day = '{date_yesterday_start}';
                        ''', engine).to_numpy()
    
    return cycles_sum,equivalent_range_mean,equivalent_range_std,crack_yesterday

def montecarlo_simulation(engine):
    
    cycles_sum,equivalent_range_mean,equivalent_range_std,crack_yesterday = load_data(engine)
    
    for j in range(len(crack_yesterday[0][0])):
        cycles_sum_rainflow_j = cycles_sum[j]
        eqst_rainflow_mean_j = equivalent_range_mean[j]
        eqst_rainflow_std_j = equivalent_range_std[j]
        crack_yesterday_mean = crack_yesterday[0][0][j]
        crack_yesterday_std = crack_yesterday[0][1][j]

        count_1 = 0;  #count failure cases
        a = np.zeros(number_of_runs)
        C = np.zeros(number_of_runs)
        da = np.zeros(number_of_runs)
        a_n = np.zeros(number_of_runs)

        for i in range(number_of_runs):
            dsig_eq_rand = np.random.normal(eqst_rainflow_mean_j,eqst_rainflow_std_j)
            a[i] = np.random.normal(crack_yesterday_mean,crack_yesterday_std) # a random generation
            C[i] = np.random.lognormal(C_mu,C_sigma) # C random generation
            da[i] = C[i]*cycles_sum_rainflow_j*(Y*dsig_eq_rand*math.sqrt(math.pi*a[i]))**m; #delta-a
            a_n[i] = a[i]+da[i]; #new crack size
            if a_n[i] >= a_c[j]:
                count_1 +=1
        P_cracks.append(count_1/number_of_runs)#Probability of failure
        a_mu_cracks.append(np.mean(a_n)) #mean new crack size
        a_sigma_cracks.append(np.std(a_n)) #standard deviation crack size
        
    df = pd.DataFrame({
                'mean_crack' : [a_mu_cracks],
                'std_crack_size' : [a_sigma_cracks],
                'pro': [P_cracks],
                'cycles_sum':[cycles_sum.tolist()],
                'monitored_day':date_yesterday_end.split()[0]
    })
    
    df.to_sql("crack_monitoring", engine, if_exists='append', index=None)

