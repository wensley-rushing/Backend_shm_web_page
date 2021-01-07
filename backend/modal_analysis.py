import datetime

import numpy as np
import pandas as pd

import backend.SSIPascoBridge.SSI_COV_AD as SSICOV

F = {}
z = {}
fs =100
ts = 6

operational_upper_bound = 10000
using_channels = 5


def get_event_data(accelerations_df, interval):
    data_event = accelerations_df[
        (accelerations_df['sampledatetime'] >= interval[0]) & (accelerations_df['sampledatetime'] <= interval[1])].copy()
    data_event.drop(columns = ['sampledatetime'], inplace = True)
    return data_event

def model_analysis(accelerations_df, interval_i,engine):
    Event_data = get_event_data(accelerations_df, interval_i)
    
    Event_data_numpy = Event_data.to_numpy()
    if Event_data_numpy.shape[0] >= fs*ts:
        
        if Event_data_numpy.shape[0] > operational_upper_bound:
            Event_data_numpy = Event_data_numpy[:operational_upper_bound,:]
            
        fn, zeta, phi, fopt, dampopt = SSICOV.SSI_COV_AD(
            Event_data_numpy[:,:using_channels],
            fs,
            ts,
            using_channels,
            30,
            5)
        
        phi_matrix = np.zeros((len(phi),using_channels))
        for i in range(len(phi)):
            phi_matrix[i][:] = np.mean(phi[i][:], axis = 1)

        df = pd.DataFrame({
            'date1' : interval_i[0],
            'date2' : interval_i[1],
            'fre' : fopt.tolist(),
            'dampopt' : dampopt.tolist(),
            'phi' : phi_matrix.tolist()
        })
        df.to_sql("modal_analysis_results", engine, if_exists='append', index=None)
        return df