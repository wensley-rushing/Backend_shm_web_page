import numpy as np
import pandas as pd
from scipy.optimize import minimize

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


last_date_query = "SELECT * FROM modal_analysis_results ORDER BY date1 DESC LIMIT 1"

#------------------------------ 7. Data base results -------------------------#
def get_model_analysis_data(engine):
    last_date = pd.read_sql(last_date_query, engine).date1[0]
    
    last_hour_query = f"SELECT * FROM modal_analysis_results WHERE date1 = '{last_date}'"
    df_model_analysis_date = pd.read_sql(last_hour_query, engine)
    
    df_model_analysis_date.sort_values(by=['fre'], ascending=False,inplace = True)
    data_numpy = df_model_analysis_date.to_numpy()
    fre = data_numpy[:,2]
    dampopt = data_numpy[:,3]
    phi = np.array(df_model_analysis_date['phi'].values.tolist()).transpose()
    return fre,dampopt,phi

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