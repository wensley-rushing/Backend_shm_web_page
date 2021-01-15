from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy import signal
import pickle


import backend.Model.Function_SHM_HAILCREEK as FHK

#------------------------------ 1.  Integration Points -----------------------#
N = 6
#------------------------------ 2.  Number of modes --------------------------#
Nm = 2
#------------------------------ 2.  Flag identifier to load cases ------------#
Ld =0
#------------------------------ 3. Semirigid Connection ----------------------#
# Semirigid Nodes
NNI = []
NNJ = []
# Semirigid Material
WW = []
Ei = []
# Flag 0 Run Zero length
flag = 0
#------------------------------ 4. Mass ---------------------------------------#
# Mass in each column simulating the mass of the hopper
Mo = 108333
# Mass in each Bin
MB1 = 0
MB2 = 0
# Mass Distribution
Mass_Case = 1
#------------------------------ 5.  Constants --------------------------------#
gama,g=78.5*10**-6,9800.
#------------------------------ 6. Intial Elastic Model ----------------------#
# Elastic modules 
E,G = 200*10**3,0.25*10**3


australian_time = datetime.now() + timedelta(hours = 10)
psd_path = Path('/home/ec2-user/PSD')

#------------------------------ 7. Data base results -------------------------#
def get_model_analysis_data(engine):
    '''
    Gets the last data from the previous day and executes the algorithms of model updating.
    
    Parameters
    ----------
    engine: sqlachemy object
        engine that connects to the database.
    
    Returns
    -------
    fre: numpy array
        array containing the frequencies of the previous day.
    dampopt: numpy array
        array containing the dampopt of the previous day.
    phi: numpy array of arrays
        Matrix containing the phi of the previous day.
    '''
    last_date_query = "SELECT * FROM modal_analysis_results ORDER BY date1 DESC LIMIT 1"
    
    last_date = pd.read_sql(last_date_query, engine).date1[0]
    last_hour_query = f"SELECT * FROM modal_analysis_results WHERE date1 = '{last_date}'"
    
    df_model_analysis_date = pd.read_sql(last_hour_query, engine)
    df_model_analysis_date.sort_values(by=['fre'], ascending=False,inplace = True)
    data_numpy = df_model_analysis_date.to_numpy()
    fre = data_numpy[:,1]
    dampopt = data_numpy[:,2]
    phi = np.array(df_model_analysis_date['phi'].values.tolist()).transpose()
    
    df_last = get_last_hour_accelerations(engine, df_model_analysis_date)
    Nfft = len(df_last[:,1]/2)+1
    df_len = int(Nfft/2+1)
    data_result = np.zeros((df_len,df_last.shape[1]+1))
    for i in range(df_last.shape[1]):
        f,pxx = signal.welch(df_last[:,i],100,nfft = Nfft,nperseg=1024)
        data_result[:,i+1] = pxx
    data_result[:,0] = f
    psd = {
    'fre':fre,
    'dampopt':dampopt,
    'df':data_result
    }
    
    file = f'{australian_time.year}-{australian_time.month:02}-{(australian_time.day-1):02}-psd.pkl'
    f = open(psd_path/file,"wb")
    pickle.dump(psd,f)
    f.close()
    return fre,dampopt,phi

def get_last_hour_accelerations(engine, df):
    '''
    Given a DataFrame object returns the accelerations of that given time.
    
    Parameters
    ----------
    engine: sqlalchemy object
        engine that connects to the database.
    df: DataFrame object
        DataFrame object containing the last event of the previous day.
    
    Returns
    -------
    accelerations_df_interval: numpy array of arrays
        Matrix containing the accelerations of the last event of the previous day.
    '''
    date1 = datetime.utcfromtimestamp(df['date1'].unique()[0].tolist()/1e9)
    date2 = datetime.utcfromtimestamp(df['date2'].unique()[0].tolist()/1e9)
    get_data_acceleration_query = f'''SELECT * FROM "100hz_python_SIM" WHERE sampledatetime BETWEEN '{date1}' AND '{date2}';'''
    acceleration_data_interval = pd.read_sql(get_data_acceleration_query, engine)
    columns = ['AI1','AI2','AI3','AI4','AI5','AI10','AI11','AI12']
    accelerations_df_interval = acceleration_data_interval[columns]
    return accelerations_df_interval.to_numpy()

#------------------------------ 4. Mass ---------------------------------------#
def Model_updating(N,Nm,E,G,Mbo,Mb1,Mb2,Mass_Case,Ld): 
    Nodes, conect, idele,Members,Names = FHK.Topology() 
    masa = FHK.Hailcreek(Members,Nodes,conect,idele, E,Mbo,Ld,N,Names,flag,NNI,WW,Ei)
    Ew = FHK.load_cases(Ld,conect,idele,Members,gama)
    dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)      
    FHK.mass_distirbution(masa,Mbo,Mb1,Mb2,dist_b1,dist_b2)

    freq,U_ops = FHK.Modal_analysis(Nm,0)
    return np.array(freq)

def opti(Freq,phi):
    
    Mbo = Mo
    
    Freq = Freq[:Nm]/30
    phi = phi[:,:Nm]
    
    MB1 = 5
    MB2 = 5
    fun = lambda x: Model_updating(N,Nm,E,G,Mbo,10**x[0],10**x[1],Mass_Case,Ld)
    def cost(Freq,x):
        Freq = np.array(Freq)
        cost = np.linalg.norm(Freq-fun(x))
        return cost
    x0 = [MB1,MB2]
    fun2 = lambda x: cost(Freq,x)
    res = minimize(fun2, x0, options = {
        'maxiter':100
    })
    return res.x

def create_pickle_visualization():
    FHK.Drawing2database(Nm)