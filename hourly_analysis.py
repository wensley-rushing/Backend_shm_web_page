import os
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

import backend.database_connection as db_con
import backend.event_detector as event_detect
import backend.time_history as time_hst
import backend.modal_analysis as modal_anl
import backend.model_updating as model_upt
import backend.rainflow as rainflow

startTime = datetime.now()


date_now = datetime.now() + timedelta(hours = 10)
year = date_now.year
month = date_now.month
day = date_now.day
hour = date_now.hour
Q2 = datetime(year, month, day, hour, 0, 0)
Q1 = Q2 - timedelta(hours = 1)

print("Don't forget Bocchi")
print(Q1)
print('Loading data')

accelerations_df, engine = db_con.connection_and_data_retrieving(Q1,Q2)

print(f'Data Loaded from Database successfully: There are {len(accelerations_df)} registries.')
print(f'Detecting Events')
events_date = event_detect.get_events(accelerations_df)

print(f'{len(events_date[1])} Events')

# Modal Analysis and Model Update
print('Doing Modal Analysis')
df_model_analysis = pd.DataFrame(columns = ['date1','date2',
                                            'fre','dampopt',
                                           'phi'])

if events_date[1].shape[0] == 1:
    df = modal_anl.model_analysis(accelerations_df, [events_date[1][0][1],
                                         np.datetime64(accelerations_df.tail(1).values[0][0])],engine)
    df_model_analysis = df_model_analysis.append(df,ignore_index=True)
else:
    for i in range(events_date[1].shape[0]-1):
        df = modal_anl.model_analysis(accelerations_df, [events_date[1][i][1],
                                         events_date[1][i+1][0]],engine)
        df_model_analysis = df_model_analysis.append(df,ignore_index=True)

    
print('Done Modal Analysis')
#time_between_events = df_model_analysis['date1'].unique()

#print('Doing Modal Update')

#for i in range(len(time_between_events)):
#    fre,dampopt,phi = model_upt.get_model_analysis_data(df_model_analysis,time_between_events[i])
#    opt = model_upt.opti(fre,phi)

#model_upt.create_pickle_visualization()
#print('Done Modal Update')

print('Doing hist analysis')
for i in range(len(events_date[1])):
    time_hst.Time_history(accelerations_df,events_date[1][i], 'AI1')
    
print('Doing rainflow')
rainflow.store_results_in_db(engine)

print('Files deleated')
engine.dispose()
print(datetime.now() - startTime)