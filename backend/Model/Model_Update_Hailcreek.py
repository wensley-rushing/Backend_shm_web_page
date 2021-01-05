#                               Hailcreek web version
#                                       2020

import numpy as np
import matplotlib.pyplot as plt
import Function_SHM_HAILCREEK as FHK
import BSDTA as  BS

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
#------------------------------ 7. Data base results -------------------------#


#------------------------------ 4. Mass ---------------------------------------#

def Model_updating(N,Nm,E,G,Mbo,Mbi,Mb2,Mass_Case,Ld):
    
    Nodes, conect, idele,Members,Names = FHK.Topology() 
    masa = FHK.Hailcreek(Members,Nodes,conect,idele, E,Mbo,Ld,N,Names,flag,NNI,WW,Ei)
    Ew = FHK.load_cases(Ld,conect,idele,Members,gama)
    dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)      
    FHK.mass_distirbution(masa,Mbo,MB1,MB2,dist_b1,dist_b2)

    freq,U_ops = FHK.Modal_analysis(Nm,0)
    return np.array(freq)

freq = Model_updating(N,Nm,E,G,Mo,MB1,MB2,Mass_Case,Ld)


def opti(Freq,phi,N,Nm,E,G,Mbo,Mass_Case,Ld):
    # breakpoint()
    MB1 = 0
    MB2 = 0
    fun = lambda x: Model_updating(N,Nm,E,G,Mbo,x[0],x[1],Mass_Case,Ld)
    def cost(Freq,x):
        # breakpoint()
        Freq = np.array(Freq)
        
        cost = np.linalg.norm(Freq-fun(x))
        return cost
    from scipy.optimize import minimize
    x0 = [MB1,MB2]
    fun2 = lambda x: cost(Freq,x)
    res = minimize(fun2, x0)
    return res.x
    
xopt = opti(freq,[],N,Nm,E,G,Mo,Mass_Case,Ld)

freq = Model_updating(N,Nm,E,G,Mo,MB1,MB2,Mass_Case,Ld)


# BS.MacVal(U[1,:].T,U[:,1].T)
# print('Mode:','Frequency [Hz]')
# for i in range(len(freq)):
#     print(i+1,round(freq[i],4))

