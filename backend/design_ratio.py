from pathlib import Path

import numpy as np

import backend.Design_Ratio.Functions_Design_ratio_HK as FHK


current_path = Path('/home/ec2-user/code_comp_files')

#------------------------------ 1.  Integration Points -----------------------#
N = 6
#------------------------------ 2.  Number of modes --------------------------#
Nm = 2
#------------------------------ 2.  Flag identifier to load cases ------------#
Ld =4
#------------------------------ 3. Semirigid Connection ----------------------#
# Semirigid Nodes
NNI = [ ]
NNJ = []
# Semirigid Material
WW = [] 
Ei = []
# Flag 0 Run Zero length
flag = 0
#------------------------------ 4. Mass --------------------------------------#
Mo = 108333
# Mass in each Bin
MB1 = 1666666.667
MB2 = 1666666.667
# Mass Distribution
Mass_Case = 1
#------------------------------ 5.  Constants --------------------------------#
gama,g=78.5*10**-6,9800.
#------------------------------ 6. Intial Elastic Model ----------------------#
# Elastic modules 
E,G = 200*10**3,0.25*10**3
# Load Factor Mass
LfM = 1.2
# Load factor for wind X
Lfwx = 0
# Load factor for wind y
Lfwy = 1


def create_file():
    Nodes, conect, idele,Members,Names = FHK.Topology()
    masa = FHK.Hailcreek(Members,Nodes,conect,idele, E,Mo,Ld,N,Names,flag,NNI,WW,Ei)
    Ew = FHK.all_loads(2,conect,Members,gama,idele,[Lfwx,Lfwy,0],LfM)
    FHK.assignloads_dist(Ew)
    # Ew = FHK.load_cases(Ld,conect,idele,Members,gama)
    dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)      
    FHK.mass_distirbution(masa,Mo,MB1,MB2 ,dist_b1,dist_b2,LfM)
    FHK.run_model()
    S,Design_ratio = FHK.Desing_Stress_Forces2(Ld,conect,idele,Ew,Names)

    code = list(Design_ratio.values())
    np.savetxt(current_path/'code.txt',code)