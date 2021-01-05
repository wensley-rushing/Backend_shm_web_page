from pathlib import Path

import numpy as np
import openseespy.postprocessing.ops_vis as opsv

import backend.TimeHistory.Functions_Hailcreek as FHK


#  integrations points
N = 6
# Number of modal frequencies
Nm = 10
# Geometry parameters
gama,g,E,G =78.5*10**-6,9800,200*10**3,0.25*10**3
# Mass in each column simulating the mass of the hopper
Mo = 108333
# Flag identifier to load cases
Ld =0
#Recorders
Mt,Tt,Et = np.array([]),np.array([]),np.array([])

Ndata = 60000

flag = 0

sw = 0
fs = 500
Mass_Case = 1

NNI = []
NNJ = []

## Change paths
Nodes, conect, idele,Members,Names = FHK.Topology() 
for k in range(len(NNI)):
    Nodes = FHK.Create_conect_copynode(NNI[k],NNJ[k],conect,Nodes)
    
def get_sensor_data(accelerations_df, interval, sensor):
    Data_sensor = accelerations_df[['sampledatetime',sensor]]
    data_interval = Data_sensor[(Data_sensor['sampledatetime'] >= interval[0]) & (Data_sensor['sampledatetime'] <= interval[1])]
    sensor_txt = sensor+'.txt'
    np.savetxt(Path.cwd()/'backend'/'data_generated'/'timehistoryfiles'/'sensors'/sensor_txt,data_interval[sensor], fmt='%s')
    
    
def Time_history(accelerations_df,interval_i, sensor):
    if sw == 0:
        WW = []
        Ei = []
        NNI = []
        Mbo = 108333 
        MB1 = 0
        MB2 = 0

    else:   
        WW = [10**20,10**20]
        Ei = [10**20,10**20]
        Mbo = 108333 
        MB1 =  9806.65*j
        MB2 =  9806.65*j 
    
    get_sensor_data(accelerations_df, interval_i, sensor)
    dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)

    #------------------------------ 8. Hailcreek----------------------------------#
    masa = FHK.Hailcreek(Members,Nodes,conect,idele, E,Mbo,Ld,N,Names,flag,NNI,WW,Ei)
    #---------------------------- 21. Mass distribution --------------------------#  
    FHK.mass_distirbution(masa,Mbo,MB1,MB2,dist_b1,dist_b2)
    #---------------------------- 9. Load Cases ----------------------------------#
    Ew = FHK.load_cases(Ld,conect,idele,Members,gama)
    #---------------------------- 14.Timeseries ----------------------------------#
    sensor_txt = sensor+'.txt'
    result_txt = sensor+np.datetime_as_string(interval_i[0]).replace(':','-')+'.txt'

    FHK.Timehistory_crusher(
        'backend/data_generated/timehistoryfiles/results/'+result_txt,
        500,
        'backend/data_generated/timehistoryfiles/sensors/'+sensor_txt)