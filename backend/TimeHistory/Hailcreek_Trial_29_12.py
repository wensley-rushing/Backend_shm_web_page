#                                   3d Hailcreek
#                                       2020

import numpy as np
import matplotlib.pyplot as plt
import Functions_Hailcreek as FHK


plt.close('all')
#                                   Global Variables

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

Simulation = 1


#Recorders
Mt,Tt,Et = np.array([]),np.array([]),np.array([])
F_M = np.zeros((Simulation,Nm))

Ndata = 60000

flag = 0

sw = 0
fs = 500
Mass_Case = 1

NNI = [ ]
NNJ = []

# NNI = [17,276]
# NNJ = [11,270]
#---------------------------- 11.Topology ------------------------------------#
Nodes, conect, idele,Members,Names = FHK.Topology() 
for k in range(len(NNI)):
    Nodes = FHK.Create_conect_copynode(NNI[k],NNJ[k],conect,Nodes)
            
#---------------------------- Cases to simulate ------------------------------#
for j in range(Simulation):
    if sw == 0:
        WW = []
        Ei = []
        NNI = []
        Mbo = 108333 
        MB1 = 0
        MB2 = 0
        dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)
     
    
    if sw == 1:   
        WW = [10**20,10**20]
        Ei = [10**20,10**20]
        Mbo = 108333 
        MB1 =  9806.65*j
        MB2 =  9806.65*j
        print('********** Mass Model *************')
        print('Mass Assign in each column M [N] ',MB1)
        #---------------------------- 22. Mass Distribution ----------------------#      

        dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)

   

        dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)      
        

             
   
    #------------------------------ 8. Hailcreek----------------------------------#
    masa = FHK.Hailcreek(Members,Nodes,conect,idele, E,Mbo,Ld,N,Names,flag,NNI,WW,Ei)
    #---------------------------- 21. Mass distribution --------------------------#  
    FHK.mass_distirbution(masa,Mbo,MB1,MB2,dist_b1,dist_b2)
    #---------------------------- 9. Load Cases ----------------------------------#
    Ew = FHK.load_cases(Ld,conect,idele,Members,gama)
    #---------------------------- 15.filecreations -------------------------------#
    # Fn = FHK.filecreations(sw,j,Mbo)
    #---------------------------- 14.Timeseries ----------------------------------#
    FHK.Timehistory_crusher('results.txt',500,'WhiteNoise.txt')



