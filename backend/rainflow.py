import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import rainflow

current_path = Path(str(__file__)).parent

path_timehistory = current_path/"data_generated"/"timehistoryfiles"/"results"


m = 3.202; #material const of the Paris law
number_of_runs = 2; #number of runs for Monte Carlo Simulation
number_analyzed_components = 8

columns = ['cycles_b1','equivalent_range_b1','cycles_b1_5','equivalent_range_b1_5',
          'cycles_b2_5','equivalent_range_b2_5','cycles_b5_5','equivalent_range_b5_5',
          'cycles_c1','equivalent_range_c1','cycles_d1','equivalent_range_d1',
          'cycles_d5','equivalent_range_d5','cycles_d1_c1','equivalent_range_d1_c1']

def crack_size(file):
    '''
    Load a file created by the time history module and executes the 
    rainflow algorithm and stores the results in an array.
    
    Parameters
    ----------
    file: txt file
        file of specific event created by the time history module.
    
    Results
    -------
    results_data: numpy array
        array with the results of the rainflow algorithm.
    '''
    stress_time_event = np.loadtxt(path_timehistory/file)
    stress_time_event_data = stress_time_event[:,:number_analyzed_components] # this most be specified at some point
    N_results = np.zeros(number_analyzed_components)
    dsig_eq_results = np.zeros(number_analyzed_components)
    results_data = np.empty((N_results.size + dsig_eq_results.size,), dtype = dsig_eq_results.dtype)
    for i in range(number_analyzed_components):
        stress_time_event_column = stress_time_event_data[:,i]
        ## RAINFLOW ALGORITHM
        n_i = list()
        dsig_i = list()
        rain_flow = rainflow.extract_cycles(stress_time_event_column); #rainflow toof the stress time history
        for range_cycle, mean, count, i_start, i_end in rain_flow:
            n_i.append(count)
            dsig_i.append(range_cycle)
        n_i = np.array(n_i);    #cycles count
        dsig_i = np.array(dsig_i); #range of stress

        ## EQUIVALENT STRESS RANGE ALGORITHM
        N_i = sum(n_i)
        dsig_eq_i = np.zeros(len(n_i))
        for j in range(len(n_i)):
            dsig_eq_i[j] = n_i[j]*dsig_i[j]**m

        dsig_eq = ((sum(dsig_eq_i))/N_i)**(1/m); #equivalent range of stress
        
        N_results[i] = N_i
        dsig_eq_results[i] = dsig_eq
    results_data[0::2] = N_results
    results_data[1::2] = dsig_eq_results
    return results_data

def store_results_in_db(engine):
    '''
    Saves the results of the crack_size functions in the database
    and then removes the files created by time history module.
    
    Parameters
    ----------
    engine: sqlachemy object
        engine that connects to the database.
    
    Returns
    -------
    None.
    '''
    files = [f for f in os.listdir(path_timehistory) if not f.startswith('.')]
    og_files = files
    results_array = np.zeros((len(files),number_analyzed_components*2))
    
    for i,file in enumerate(files):
        results_array[i] = crack_size(file)
        
    files = [file.replace('.txt', '') for file in files]
    files = [file.replace('AI1', '') for file in files]
    
    files_datetime = list()
    for file in files:
        file = file.split('T')
        text = file[1].replace("-",":")
        date_time_str = file[0]+" "+text[:15]
        files_datetime.append(datetime.strptime(date_time_str, '%Y-%m-%d %H:%M:%S.%f'))
        
    df = pd.DataFrame(results_array, columns = columns)
    df['timestamp'] = files_datetime
    df.to_sql("stress_events_rainflow", engine, if_exists='append', index=None)
    remove_files(og_files)

def remove_files(files):
    '''
    Removes the files created in the time history module because
    they are not longer needed.
    
    Parameters
    ----------
    files: list of files
        files created in the time history module.
        
    Returns
    -------
    None.
    '''
    for file in files:
        path_to_file = path_timehistory/file
        os.remove(path_to_file)