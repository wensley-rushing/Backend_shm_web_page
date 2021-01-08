import numpy as np 
import Functions_Design_ratio_HK as FHK
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
LfM = 1.25
# Load factor for wind X
Lfwx = 0
# Load factor for wind y
Lfwy = 1
combo = [[1.25,0,1,0,0,0,0],
          [1.25,0,-1,0,0,0,0],
          [1.25,1,0,0,0,0,0],
          [1.25,-1,0,0,0,0,0],
          [1.35,0,0,0,0,0,0],
          [0.5,0,0,1,0.3,0,0],
          [0.5,0,0,-1,-0.3,0,0],
          [0.5,0,0,1,-0.3,0,0],
          [0.5,0,0,-1,0.3,0,0],
          [0.5,0,0,0.3,1,0,0],
          [0.5,0,0,-0.3,-1,0,0],
          [0.5,0,0,0.3,-1,0,0],
          [0.5,0,0,0.3,1,0,0],
          [1.2,0,0,0,0,0,1.5]]
results = {}
#------------------------------ 7. Functions ---------------------------------#
Nodes, conect, idele,Members,Names = FHK.Topology()
for i in range(len(combo)):

    LfM,Lfwx,Lfwt = combo[i][:3]
    Lsx,Lsy,Lsz = combo[i][3],combo[i][4],combo[i][5]
    lT = combo[i][6]
    masa = FHK.Hailcreek(Members,Nodes,conect,idele, E,Mo,Ld,N,Names,flag,NNI,WW,Ei)
    Ew = FHK.all_loads(2,conect,Members,gama,idele,[Lfwx,Lfwy,0],LfM)
    FHK.assign_point_seismic([Lsx,Lsy,Lsz])
    FHK.assign_point_truck(lT)
    FHK.assignloads_dist(Ew)
    # Ew = FHK.load_cases(Ld,conect,idele,Members,gama)
    dist_b1,dist_b2 = FHK.mass_shapedistribution(Mass_Case)      
    FHK.mass_distirbution(masa,Mo,MB1,MB2 ,dist_b1,dist_b2,LfM)
    FHK.run_model()
    S,Design_ratio = FHK. Desing_Stress_Forces2(Ld,conect,idele,Ew,Names) 
    code = list(Design_ratio.values())
    results[i] = code
#-------------------  8. Maximun values of Load Cases Dr ---------------------#
Results = np.array(list(results.values())).T   
Results_Wind =np.max(Results[:,:3],axis = 1)
Results_Op = Results[:,4]
Results_Seismic =np.max(Results[:,-7:],axis = 1)
Results_Truck =Results[:,-1]
Max_Design_Ratio = np.vstack((Results_Wind,Results_Op,Results_Seismic,Results_Truck)).T
cases = [1,2,3,4]
lk = np.zeros((Max_Design_Ratio.shape[0],4))
for i in range(len(cases)):
    prob = FHK.Realiability(Max_Design_Ratio[:,i])
    lko = FHK.Likelyhoood(prob,cases[i])
    lk[:,i] = lko[:,0]
#------ Geting Max values ----#
lk  = np.max(lk,axis = 1)

#------ MaxDesing Ratio ----#
Code = Results_Op = Results[:,4]
np.savetxt('code.txt',Code)  
idd = np.array([23,18,6,69,184,236,262])-1
like = lk[idd]
np.savetxt('Likelihood.txt',like)