import numpy as np
import matplotlib.pyplot as plt
import opensees.openseespy as ops
import openseespy.postprocessing.ops_vis as opsv
import backend.Model.BSDTA as BS
from scipy import stats
from pathlib import Path

current_path = Path(str(__file__)).parent

def Topology():
    # Get the Nodes coordinate
    Nodes=np.genfromtxt(current_path/"Nodes.txt",invalid_raise=False)
    # Member connections 
    conect = np.loadtxt(current_path/"Connc.txt",delimiter=',')
    # Frames asigns
    idele=np.genfromtxt(current_path/"Elementid.txt",invalid_raise=False)
    # Members Properties
    Members=np.genfromtxt(current_path/"Members.txt",invalid_raise=False)
    # Names
    Names = np.genfromtxt(current_path/"Nmes.txt",dtype='str',delimiter = '')
    return Nodes, conect, idele,Members,Names

#---------------------------- 1. Model forces --------------------------------#
def Hailcreek(Members,Nodes,conect,idele, E,Mb,Ld,N,Names,flag,NNI,WW,Ei):
    ops.wipe()
    # Gmod
    G = 0.25*10**3
    # Gravity 
    g = 9800
    # Gama
    gama = 78.5*10**-6
    # Model DOf
    ops.model('basic','-ndm',3,'-ndf',6)
    # Name id sections
   
    # Nodes with constrain
    Supp = [1,4,9,14,29,34,65,67,89,94,122,125,130,135,139,142,147\
            ,152,170,175,207,209,230,235,260,263,268,273]

    # Z global axis reference
    z = [0,0,1]
    
    # Alocate memory for mass 
    masa = np.zeros((len(Nodes),2))
    # Asign each row with a Node
    for i in range(len(Nodes)):
        ops.node(i+1,Nodes[i,0],Nodes[i,1],Nodes[i,2])
        masa[i, 0] = i+1
    # 
    for i in range(len(conect[:,1])):
        
        # Node i&j Coordinate
        XYZI = ops.nodeCoord(int(conect[i,0]))
        XYZJ = ops.nodeCoord(int(conect[i,1]))
        # xaxis is parales to the element axis
        xaxis = np.subtract(XYZJ,XYZI)  
        vecxz = np.cross(xaxis,z)
        # Length of elements
        L =np.linalg.norm(xaxis)
        # Mass for each sections
        m =Members[int(idele[i])-1,3]*L*gama/2/g
       
    
        if np.linalg.norm(vecxz) == 0:
    
            ops.geomTransf('Linear', i,0,-1,0)
            Internal_Forces(i,E,G,Members,idele,N) 
            con = [int(conect[i,0]),int(conect[i,1])]       
            ops.element('forceBeamColumn',i+1,*con,i,i+1)
            index = np.where(masa[:,0] == conect[i,0])[0][0]
            masa[index,1] += m
            index = np.where(masa[:,0] == conect[i,1])[0][0]
            masa[index,1] += m
        else:
            
            if XYZI[2]-XYZI[2] == 0:
      
                ops.geomTransf('Linear', i,*vecxz)
                Internal_Forces(i,E,G,Members,idele,N)
                con = [int(conect[i,0]),int(conect[i,1])]       
                ops.element('forceBeamColumn',i+1,*con,i,i+1)
                index = np.where(masa[:,0] == conect[i,0])[0][0]
                masa[index,1] += m
                index = np.where(masa[:,0] == conect[i,1])[0][0]
                masa[index,1] += m
            else:
                Internal_Forces(i,E,G,Members,idele,N)
                ops.geomTransf('Linear', i,*vecxz)
                Internal_Forces(i,E,G,Members,idele.N)
                con = [int(conect[i,0]),int(conect[i,1])]       
                ops.element('Truss',i+1,*con,i,i+1)
                index = np.where(masa[:,0] == conect[i,0])[0][0]
                masa[index,1] += m
                index = np.where(masa[:,0] == conect[i,1])[0][0]
                masa[index,1] += m
                
                
        # print('Section  ',Names[int(idele[i])-1],'  assigned to element No.  ',i+1)
        
    # Nodes whit supports
    for i in range(len(Supp)):
        ops.fix(int(Supp[i]),1,1,1,0,0,0)
    
    for i in range(len(masa)):
        ops.mass(masa[i, 0], masa[i, 1], masa[i, 1], masa[i, 1], 0, 0, 0)
     
    if flag == 0:
        for i in range(len(NNI)):   
            zero_length(NNI[i],277+i,WW[i],Ei[i])    
                
    
    """
    In this part of the code, the constrain is created 
    
    """
    Const2 = [155, 212, 276, 271, 266, 238, 178, 145, 150]
    Const1 = [17,70,138,133,128,97,37,7,12]
   
    Cord_N12  = ops.nodeCoord(12)
    Cord_N17 = ops.nodeCoord(17)
    Cord_N138 = ops.nodeCoord(138)
    
    L_17_138 = np.subtract(Cord_N138,Cord_N17)
    L_medios = np.linalg.norm(L_17_138)/2
    Cord_Maestro = Cord_N12 
    Cord_Maestro[0] = Cord_Maestro[0]+L_medios 
    ops.node(1776,*Cord_Maestro)
    ops.fix(1776,0,0,1,1,1,0)
    ops.rigidDiaphragm(3,1776, *Const1)
    
    
    
    
    Cord_N150  = ops.nodeCoord(150)
    Cord_N155 = ops.nodeCoord(155)
    Cord_N276 = ops.nodeCoord(276)
    L_155_276 = np.subtract(Cord_N155,Cord_N276)
    L_medios = np.linalg.norm(L_155_276)*0.5
    Cord_Maestro = Cord_N150
    Cord_Maestro[0] = Cord_Maestro[0]+L_medios 
    ops.node(1777,*Cord_Maestro)
    ops.fix(1777,0,0,1,1,1,0)
    ops.rigidDiaphragm(3,1777, *Const2)
    
    ops.constraints('Transformation')
    return masa
#---------------------------- 2. Integration Points --------------------------#
def Internal_Forces(i,E,G_mod,Members,idele,N):
    #This function creates: 
    #the recorders for each fucntions
    #the elastic section 
    #section('Elastic', secTag, E_mod, A, Iz, Iy, G_mod, Jxx)
    #the integrations points with Lobatos
    #beamIntegration('Lobatto', tag, secTag, N)
    
    A = Members[int(idele[i])-1,3]
    Iz = Members[int(idele[i])-1,1]
    Iy = Members[int(idele[i])-1,2]
    Jxx = Members[int(idele[i])-1,0]
    
    ops.section('Elastic',i+1, E, A, Iz, Iy, G_mod, Jxx)
    ops.beamIntegration('Lobatto', i+1, i+1, N)  
#---------------------------- 3. Add Mass columns ----------------------------#
def Mass_bin(Const1,Const2,masa,Mb):
    """
    Mb: Point Mass in each column modeling
    the mass of the upper part of the beam
    luis value ; 108333N
    """
    g = 9800
    for i in range(len(Const1)):
        jj = int(Const1[i])-1
        M  = masa[jj,1]+Mb/g
        ops.mass(masa[jj,0],M,M,M,0,0,0)

    for i in range(len(Const2)):
        jj = int(Const2[i])-1
        M  = masa[jj,1]+Mb/g
        ops.mass(masa[jj,0],M,M,M,0,0,0)          
#---------------------------- 4. Modal Analysis   ----------------------------#       
def Modal_analysis(Nm,sw):
    """
    Nm: Number of Modes that the code is going to return 
    sw: Switch value for returning coordinates 
    """
    Const2 = [155, 212, 276, 271, 266, 238, 178, 145, 150]
    Const1 = [17,70,138,133,128,97,37,7,12]
    freq = ops.eigen(Nm)
    Nodes_AC= np.concatenate((Const1,Const2))
    
    U_ops = np.zeros((36,Nm))
    for i in range(0,Nm):
        for j in range(0,len(Nodes_AC)*2):
            if j<18:
                U_ops[j,i]= ops.nodeEigenvector(int(Nodes_AC[j]),i+1)[0]
            else:
                U_ops[j,i]= ops.nodeEigenvector(int(Nodes_AC[j-18]),i+1)[1]
    
    
    for i in range(0,Nm):
         U_ops[:,i] =  U_ops[:,i]/np.linalg.norm(U_ops[:,i])
    
    # print('Mode:','Frequency [Hz]')
    for i in range(Nm):
    
       freq[i]=freq[i]**0.5/(2*np.pi)
    
    if sw == 1:
        CORR = np.zeros((len(Nodes_AC),4))
        for i  in range(len(Nodes_AC)):
            CORR[i,0] = Nodes_AC[i]
            CORR[i,1:5] = ops.nodeCoord(int(Nodes_AC[i]))
            print(CORR[i,:])
            return freq,U_ops,CORR
    else:
        CORR =[]
        return freq,U_ops

#---------------------------- 5. Load Cases ----------------------------------#
def load_cases(ld,conect,idele,Members,gama): 
    
    """
    in this section the loads for each case is defined ,
    To do : we need to ad the wind and seismic case
    To do : remove the for loop for the mass, stored 
            the information in a variable in the 
            first asignation 
    To do: this could work with a      or with a dictionary 
    
    Some documentation 
    eleLoad('-ele', *eleTags, '-range', eleTag1, eleTag2, 
            '-type', '-beamUniform',Wy, Wz=0.0, Wx=0.0, '-beamPoint', 
            Py, Pz=0.0, xL, Px=0.0, '-beamThermal',*tempPts)    
    """
    if ld == 3:
        ops.timeSeries('Linear',1)
        ops.pattern("Plain", 1, 1)
        Ew= {}
        loads = np.loadtxt('Loads_wx.txt')
        xo = [1,0,0]
        # Y = [1,0,0]
        # Z = [0,1,0]
        for i in range(len(loads[:,0])):
            tag = loads[i,0]
            yi = ops.eleResponse(tag, 'ylocal')
            zi = ops.eleResponse(tag, 'zlocal')
            if np.linalg.norm(np.cross(xo,yi)) == 0:
                    ter = np.sum(xo)/np.sum(yi)
                    wi = loads[i,1]*ter
                    ops.eleLoad('-ele',tag,'-type','-beamUniform', wi,0,0)
                    Ew[tag] = {tag: ['-beamUniform',wi,0,0]}
            else:
                    print(tag)
                    tag = loads[i,0]
                    ter = np.sum(xo)/np.sum(zi)
                    wi = loads[i,1]*ter
                    print('¨*_****',wi,'****')
                    ops.eleLoad('-ele',tag,'-type','-beamUniform', 0,wi,0)
                    Ew[tag] = {tag: ['-beamUniform',0, wi,0]}
        """ System Properties """
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.algorithm('NewtonRhapson')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1) 
        return Ew
    
    elif  ld == 2:
        """
        In this part we put the mass of the element as load ,
        for columns the load is applies along the local x axis
        and for beams o braces the load is applied along the
        local y axis
        we need to check this !!!!!!!!!!!!!!!!! alejandro
        """
        ops.timeSeries('Linear',1)
        ops.pattern("Plain", 1, 1)
        Ew = {}
        z = [0,0,1]
        for i in range(len(conect[:,1])):
            # Node i&j Coordinate
            XYZI = ops.nodeCoord(int(conect[i,0]))
            XYZJ = ops.nodeCoord(int(conect[i,1]))
            # xaxis is parales to the element axis
            xaxis = np.subtract(XYZJ,XYZI)  
            vecxz = np.cross(xaxis,z)
            # Node i&j Coordinate
            # xaxis is parales to the element axis
            # Length of elements
            # Mass for each sections
            m =Members[int(idele[i])-1,3]*gama
            if np.linalg.norm(vecxz) == 0:
                ops.eleLoad('-ele',i+1,'-type','-beamUniform',0, m,0) 
                
                Ew[i+1] = {i+1: ['-beamUniform',0, 0,m]}
                
            else:
                ops.eleLoad('-ele',i+1,'-type','-beamUniform',m,0,0)  
                Ew[i+1] = {i+1: ['-beamUniform',m, 0,0]}
       
        """ System Properties """
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.algorithm('NewtonRhapson')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1) 
        return Ew
    if ld == 4:
        ops.timeSeries('Linear',1)
        ops.pattern("Plain", 1, 1)
        Ew= {}
        loads = np.loadtxt('Loads_wy.txt')
        xo = [0,1,0]
        # Y = [1,0,0]
        # Z = [0,1,0]
        for i in range(len(loads[:,0])):
            tag = loads[i,0]
            yi = ops.eleResponse(tag, 'ylocal')
            zi = ops.eleResponse(tag, 'zlocal')
            if np.linalg.norm(np.cross(xo,yi)) == 0:
                    ter = np.sum(xo)/np.sum(yi)
                    wi = loads[i,1]*ter
                    ops.eleLoad('-ele',tag,'-type','-beamUniform', wi,0,0)
                    Ew[tag] = {tag: ['-beamUniform',wi,0,0]}
            else:
                    print(tag)
                    tag = loads[i,0]
                    ter = np.sum(xo)/np.sum(zi)
                    wi = loads[i,1]*ter
                    print('¨*_****',wi,'****')
                    ops.eleLoad('-ele',tag,'-type','-beamUniform', 0,wi,0)
                    Ew[tag] = {tag: ['-beamUniform',0, wi,0]}
        """ System Properties """
        ops.system('BandSPD')
        ops.numberer('RCM')
        ops.algorithm('NewtonRhapson')
        ops.integrator('LoadControl', 1)
        ops.analysis('Static')
        ops.analyze(1) 
        return Ew
    if ld == 5:
        ex = '.txt'
        root2 = 'Global_internal_forces/'
        ops.recorder('Element', '-file',root2 +str(323)+ex,'-ele', 323,'localForce')
            
        ops.recorder('PVD', 'Dis_pvd', '-dT', 0.01, 'disp')
        ops.timeSeries('Path', 1, '-dt', 1/50, '-filePath', 'wn.txt', '-factor',1)
        ops.pattern('Plain', 1, 1, '-accel', 1)    
        ops.sp(70, 1, 1)
        ops.sp(212, 1, 1)
        #Rayleigh
        freq = ops.eigen(2)
        dampRatio = 0.0021
        ops.rayleigh(0., 0., 0., 2*dampRatio/freq[0])
        # create the analysis
   
        ops.numberer('Plain')    # renumber dof's to minimize band-width (optimization), if you want to
        ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
        ops.algorithm('Newton')     # use Linear algorithm for linear analysis
        ops.integrator('Newmark', 0.5, 0.25)    # determine the next time step for an analysis
        ops.analysis('Transient')   # define type of analysis: time-dependent
        ops.analyze(1000, 1/50)     # apply 3995 0.01-sec time steps in analysis
         
        ops.remove('recorders')
        Ew = {0}
    if ld == 7:
        ex = '.txt'
        root2 = 'Global_internal_forces_800T_100t/'
        ops.recorder('Element', '-file',root2 +str(6)+ex,'-ele',12,'localForce')
        Ew ={}
        ops.recorder('PVD', 'acc_pvd', '-dT', 1/7, 'disp')
        ops.timeSeries('Path', 1, '-dt', 1/14, '-filePath', 'Loads_800T.txt', '-factor',1)
        ops.pattern('Plain', 1,1)    
        ops.load(1776,1,0,0,0,0,0)
        freq = ops.eigen(2)
        dampRatio = 0.0021
        ops.rayleigh(0., 0., 0., 2*dampRatio/freq[0])
        ops.numberer('Plain')    # renumber dof's to minimize band-width (optimization), if you want to
        ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
        ops.algorithm('Linear')     # use Linear algorithm for linear analysis
        ops.integrator('Newmark', 0.5, 0.25)    # determine the next time step for an analysis
        ops.analysis('Transient')   # define type of analysis: time-dependent
        ops.analyze(250, 1/7)     # apply 3995 0.01-sec time steps in analysis
        ops.remove('recorders')
        return Ew
    if ld  == 0:
        Ew= {}
        return Ew
#---------------------------- 6. Create hinges ------------------------------# 

def Create_conect_copynode(Ni,Nj,conect,Nodes):
    
    """
    Ni = node to copy
    Nj = node to connect the copy
    Connect = conection array
    Nodes = coordinates
    """

    a = np.where(conect[:,0]== Ni)[0][:]
    b = np.where(conect[:,1]== Ni)[0][:]

    sw = 0
    while sw == 0:
        for k in range(len(a)):
            if conect[a[k],1] == Nj:
                sw =1 
    
                # print('creo copia de i y uno a Nj', '*¨ ',conect[a[k]])
                Nodes = np.vstack((Nodes,Nodes[Ni-1,:]))
                conect[a[k],0] = int(len(Nodes[:,1]))
                break
            else:
                sw = 0
        if sw == 0:
            for k in range(len(b)):
                if conect[b[k],0] == Nj:
                    sw =1 
                    Nodes = np.vstack((Nodes,Nodes[Ni-1,:]))
                    conect[b[k],1] = int(len(Nodes[:,1]))
                    break
                else:
                    sw = 0
                    print('Couldnt find the node')
                
    return Nodes
#---------------------------- 7.Topology ------------------------------------#
def zero_length(Ni,Nic,EO,Ei):
    taag =int( Nic+1000)
    eleNodes = [Nic,Ni]
    N1,N2 = int(taag+10),int(taag*9)
    ops.uniaxialMaterial('Elastic',int(taag+10),Ei)
    ops.uniaxialMaterial('Elastic',int(taag*9),Ei)
    ops.uniaxialMaterial('Elastic',taag, EO)
    Materials = [N1,N2,N1,N1,N2,N1]
    #ops.equalDOF(eleNodes[1], eleNodes[0], *[4,5,6])
    ops.element('zeroLength',taag , *eleNodes, '-mat', *Materials, '-dir', *[1,2,3,4,5,6])
    
#---------------------------- 8. Mass d4.istribution --------------------------#  
def mass_distirbution(masa,Mbo,MB1,MB2,dist_b1,dist_b2):

    MB1 = 9*MB1
    MB2 = 9*MB2
     
    Const2 = [155, 212, 276, 271, 266, 238, 178, 145, 150]
    Const1 = [17,70,138,133,128,97,37,7,12]
    """
    Mb: Point Mass in each column modeling
    the mass of the upper part of the beam
    luis value ; 108333N
    """
    g = 9800
    for i in range(len(Const1)):
        Mb = Mbo+MB1*dist_b1[i]
        jj = int(Const1[i])-1
        M  = masa[jj,1]+Mb/g
        ops.mass(masa[jj,0],M,M,M,0,0,0)

    for i in range(len(Const2)):
        Mb = Mbo+MB2*dist_b2[i]
        jj = int(Const2[i])-1
        M = masa[jj,1]+Mb/g
        ops.mass(masa[jj,0],M,M,M,0,0,0)   

#---------------------------- 22. Mass Distribution --------------------------#     


def mass_shapedistribution(Mass_Case):
    if Mass_Case == 1:
        x1 = 1/9
        dist_b1 = [x1,x1, x1, x1, x1,x1, x1, x1,x1]
        dist_b2 = [x1,x1, x1, x1, x1,x1, x1, x1,x1]
        return dist_b1,dist_b2
    if Mass_Case == 2:
 
        x1 = 0.2
        x2 = 0.3
        x3 = 1+(3*0.8+2*0.5)/4
        dist_b1 = [x1,x1, x1, x2, x3,x3, x3, x3,x2]
        dist_b2 = [x1,x1, x1, x2, x3,x3, x3, x3,x2]
        return dist_b1,dist_b2
    if Mass_Case == 3:
        x1 = 1.1
        x2 = 1
        dist_b1 = [x1,x2, x1,x2, x1,x2, x2, x1,x2]
        dist_b2 = [x1,x2, x1,x2, x1,x2, x2, x1,x2]
        return dist_b1,dist_b2

#------------------------------- 35. Get Modeshapes --------------------------#    

def getmodeshapes(Nodes,Nm):
    Nd = Nodes.shape[0]
    Modes = np.zeros([Nd,3])
    for i in range(Nd):
        Modes[i,0]=ops.nodeEigenvector(int(i+1), Nm)[0]
        Modes[i,1]=ops.nodeEigenvector(int(i+1), Nm)[1]
        Modes[i,2]=ops.nodeEigenvector(int(i+1), Nm)[2]
        
    return Modes

def create_pickles(Nodes,conect,path,name):
    import pickle
    Xe = []#np.zeros((len(conect[:,0]),2))
    Ye = []#np.zeros((len(conect[:,0]),2))
    Ze = []#np.zeros((len(conect[:,0]),2))
    for i in range(len(conect[:,0])):
        ii = conect[i,0]
        jj = conect[i,1]
        corrii = Nodes[int(ii-1),:]
        corrjj = Nodes[int(jj-1),:]
        Xe += [corrii[0],corrjj[0],'None']
        Ye += [corrii[1],corrjj[1],'None']
        Ze += [corrii[2],corrjj[2],'None']
    kex = name+"Xe.pkl"
    key = name+"Ye.pkl"
    kez = name+"Ze.pkl"
    f = open(path/kex,"wb")
    pickle.dump(Xe,f)
    f.close()
    f = open(path/key,"wb")
    pickle.dump(Ye,f)
    f.close()
    f = open(path/kez,"wb")
    pickle.dump(Ze,f)
    f.close()

def Drawing2database(Nm):
    ## Ask Alejandro why is always returnig the same
    FS = 50000
    Nodes, conect, idele,Members,Names = Topology()
    path = current_path/".."/".."/".."/"pickle_visualization"
    for i in range(Nm):
        idd  = str(i+1)
        Modes = getmodeshapes(Nodes,i+1)
        Nodes2 = Nodes+Modes*FS
        create_pickles(Nodes2,conect,path,f"Mode{idd}")
        name = f"Mode{idd}.txt"
        np.savetxt(path/f"Mode{idd}.txt",Nodes2)