#                                   3d Hailcreek
#                                       2020

import numpy as np
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import openseespy.postprocessing.ops_vis as opsv
# import BSDTA as BS
from scipy import stats
from pathlib import Path

current_path = Path(str(__file__)).parent


"""
    1. rot_transf_3d
        Coordinate Transformations
    2. Section_forces_distributio_3d
        Get internal forces
    3. Internal Foces: 
        Create internal forces
    4. Mass_bin
        Add Mass to the columns
    5. Modal_analysis
        return Frequencies and modeshapes
    6. Forces 
        retrieve teh internal forces
        each element
    7. Desing_Stress_Forces
        return de stress of each element
    8. Hailcreek
        Create de opensees model 
    
        
-------------  $$ ------------------------   
     Ld: Load case identifier
    1. Ld = 2
        Convert Mass to a load cases
    2. Ld = 3
        Wind Load acting X axis
    3. Ld = 4
        Wind Load acting Y axis
    4. Ld = 5 
        Time series analysis, acceleration 
        record acting in the middle
        of the constrains
    5. Ld = 6
        Time series analysis, Forces
        acting in the middle of
        the constrains 
        
"""


#---------------------------- 1.Rot Tranf ------------------------------------#
def rot_transf_3d(ex, ey, ez, g):

    Lxyz = np.array([ex[1]-ex[0], ey[1]-ey[0], ez[1]-ez[0]])

    L = np.sqrt(Lxyz @ Lxyz)

    z = np.zeros((3, 3))

    G = np.block([[g, z, z, z], 
                  [z, g, z, z],
                  [z, z, g, z],
                  [z, z, z, g]])

    return G, L
#---------------------------- 2. Internal Forces------------------------------#
def section_force_distribution_3d(ex, ey, ez, pl, nep=17,
                                  ele_load_data=['-beamUniform', 0., 0., 0.]):
    """
    Calculate section forces (N, Vy, Vz, T, My, Mz) for an elastic 3d beam.

    Longer description

    Parameters
    ----------

    ex : list
        x element coordinates
    ey : list
        y element coordinates
    ez : list
        z element coordinates
    pl : ndarray
    nep : int
        number of evaluation points, by default (2) at element ends

    ele_load_list : list
        list of transverse and longitudinal element load
        syntax: [ele_load_type, Wy, Wz, Wx]
        For now only '-beamUniform' element load type is acceptable.

    Returns
    -------

    s : ndarray
        [N Vx Vy T My Mz]; shape: (nep,6)
        column vectors of section forces along local x-axis

    uvwfi : ndarray
        [u v w fi]; shape (nep,4)
        displacements at nep points along local x

    xl : ndarray
        coordinates of local x-axis; shape (nep,)

    Notes
    -----

    Todo: add '-beamPoint' element load type

    """

    # eload_type = ele_load_data[0]
    Wy, Wz, Wx = ele_load_data[1], ele_load_data[2], ele_load_data[3]

    N1, Vy1, Vz1, T1, My1, Mz1 = pl[:6]

    Lxyz = np.array([ex[1]-ex[0], ey[1]-ey[0], ez[1]-ez[0]])
    L = np.sqrt(Lxyz @ Lxyz)

    xl = np.linspace(0., L, nep)
    one = np.ones(nep)

    N = -1.*(N1*one + Wx*xl)
    Vy = Vy1*one + Wy*xl
    Vz = Vz1*one + Wz*xl
    T = -T1*one
    Mz = -Mz1*one + Vy1*xl + 0.5*Wy*xl**2
    My = My1*one + Vz1*xl + 0.5*Wz*xl**2

    s = np.column_stack((N, Vy, Vz, T, My, Mz))

    return s, xl
#---------------------------- 3. Integration Points --------------------------#
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
    
#---------------------------- 4. Add Mass columns ----------------------------#
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
        
#---------------------------- 5. Modal Analysis ..----------------------------#
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
    
    print('Mode:','Frequency [Hz]')
    for i in range(Nm):
        print(i+1,freq[i]**0.5/(2*np.pi))
    
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

#---------------------------- 6. Moment Diagrams -----------------------------#
def Forces(ele_tag,Ew):
    
    """
    here in force i need to run an analysis 
    for each load case for is only implemented
    the load cases that uses the mass as a load    
    
    """
    
    
    nd1, nd2 = ops.eleNodes(ele_tag)
    # element x, y coordinates
    ex = np.array([ops.nodeCoord(nd1)[0],
                   ops.nodeCoord(nd2)[0]])
    ey = np.array([ops.nodeCoord(nd1)[1],
                   ops.nodeCoord(nd2)[1]])
    ez = np.array([ops.nodeCoord(nd1)[2],
                   ops.nodeCoord(nd2)[2]])
    
    # eo = Eo[i, :]
    xloc = ops.eleResponse(ele_tag, 'xlocal')
    yloc = ops.eleResponse(ele_tag, 'ylocal')
    zloc = ops.eleResponse(ele_tag, 'zlocal')
    g = np.vstack((xloc, yloc, zloc))
    G, _ = rot_transf_3d(ex, ey, ez, g)
    g = G[:3, :3]
        
    g = np.vstack((xloc, yloc, zloc))
    pl = ops.eleResponse(ele_tag, 'localForces')
    s,x = section_force_distribution_3d(ex, ey, ez, pl,17,Ew[ele_tag])
    Mz = []
    Mz= np.append(Mz,0)
    Mz = np.append(Mz,s[:,5])
    Mz= np.append(Mz,0)
    
    
    return s,Mz

#---------------------------- 7. Desing Stress Forces ------------------------#
def Desing_Stress_Forces(ld,conect,idele,Ew):
    Ful_Sec = np.loadtxt('Members_full.txt')
    if ld<2:
        for i in range(len(conect[:,1])):
            #(N, Vy, Vz, T, My, Mz)
            tag = i+1
            s,Mz =Forces(tag,Ew[tag])
    
            print('Element id:',tag,'max Vy: ',max(s[:,1]))
            print('Element id:',tag,'min Vy: ',min(s[:,1]))
            print('Element id:',tag,'max P: ',max(s[:,0]))
            print('Element id:',tag,'min P: ',min(s[:,0]))
            
            if tag == 2:
                np.savetxt('2.txt',s,delimiter = '  ')
    else:       
        if ld >= 5:
            return 
        else:
            if ld == 3:
                Loads = np.loadtxt('Loads_wx.txt')
            else:
                Loads = np.loadtxt('Loads_wy.txt') 
                      
            for i in range(len(Loads[:,1])):
                #(N, Vy, Vz, T, My, Mz)
                tag = Loads[i,0]
                s,Mz =Forces(tag,Ew[tag])
                Zy = Ful_Sec[int(idele[int(tag-1)]-1),15]
                Zx = Ful_Sec[int(idele[int(tag-1)]-1),16]
                AA = Ful_Sec[int(idele[int(tag-1)]-1),6]
                Asy = Ful_Sec[int(idele[int(tag-1)]-1),12]
                Asx = Ful_Sec[int(idele[int(tag-1)]-1),11]
                print( 'Zy',' Section: ' , tag ,': ',Zy)
                print( 'Zx',' Section: ' , tag ,': ',Zx)
                print('Elemento id:',tag,'max Mz: ',max(Mz)/1000000)
                print('Elemento id:',tag,'min Mz: ',min(Mz)/1000000)
                print('Element id:',tag,'max My: ',max(s[:,4])/1000000)
                print('Element id:',tag,'min My: ',min(s[:,4])/1000000)
                print('****** Stress *******')
                print('Element id:',tag,'max Sigma My + : ',max(s[:,5])/Zy)
                print('Element id:',tag,'min Sigma My - : ',min(s[:,5])/Zy)
                print('Element id:',tag,'max Sigma My + : ',max(s[:,4])/Zx)
                print('Element id:',tag,'min Sigma My - : ',min(s[:,4])/Zx)
                print('Element id:',tag,'Normal Stress  : ',max(s[:,0])/AA)
                print('Element id:',tag,'Tau y +  : ',max(s[:,2])/Asy)
                print('Element id:',tag,'Tau y -  : ',min(s[:,2])/Asy)
                print('Element id:',tag,'Tau x +  : ',max(s[:,1])/Asx)
                print('Element id:',tag,'Tau x -  : ',min(s[:,1])/Asx)
                
                print('****** Forces *******')
                
                if tag == 496:
                    np.savetxt('496.txt',s,delimiter = '  ')
                    
#---------------------------- 8. Model forces --------------------------------#
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

#---------------------------- 9. Load Cases ----------------------------------#
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
    
#---------------------------- 10. Create hinges ------------------------------# 
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

#---------------------------- 11.Topology ------------------------------------#
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
#---------------------------- 12.Topology ------------------------------------#
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
    
#---------------------------- 13.Temperature Model ---------------------------#
def Temperaturemodel(Tf):
        To = 20
        dt = Tf-To
        e1 = 3.768
        e2 = 1.0
        e3 = 639.0
        e4 = 1650
        Eo =  200*10**3
        E = Eo*(np.exp(-0.5*(dt/e3)**e1-0.5*(dt/e4)**e2))
        return E
    
#---------------------------- 14.Timeseries ----------------------------------#
def Timehistory(Fn,Ndata,name):
    
    Nrsp = [155, 212, 276, 271, 266, 238, 178, 145, 150,  17,  70, 138, 133,\
       128,  97,  37,   7,  12]
    ops.wipeAnalysis()
    freq = ops.eigen(4)
    ops.timeSeries('Path', 1, '-dt', 1/100, '-filePath', name, '-factor',1)
    ops.recorder('Node', '-file', Fn  ,'-timeSeries', 1,'-dT', 1/250,'-node', *Nrsp,'-dof',1,'accel')
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)    
    dampRatio = 0.002
    ops.rayleigh(0., 0., 0., 2*dampRatio/freq[0])
    # create the analysis
    ops.constraints(' Transformation')
    ops.numberer('Plain')    # renumber dof's to minimize band-width (optimization), if you want to
    ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
    ops.algorithm('Linear')     # use Linear algorithm for linear analysis
    ops.integrator('Newmark', 0.5, 0.25)    # determine the next time step for an analysis
    ops.analysis('Transient')   # define type of analysis: time-dependent
    ops.analyze(int(Ndata), 1/50)     # apply time steps in analysis
    ops.remove('recorders')

#---------------------------- 15.filecreations -------------------------------#
def filecreations(sw,i,Mss):
    if sw == 1:
        Nme ='test_'
        Mss = 'Mass'
        root = 'Acc_record/'
        ext = '.txt'
        filen = Nme+str(i)+Mss
        Fn = root+filen+ext
        return Fn
    if sw==2:
        Nme ='test_'
        Mss = 'Mass'
        root = 'Acc_record2/'
        ext = '.txt'
        filen = Nme+str(i)+Mss
        Fn = root+filen+ext
        return Fn      
    if sw==3:
        Nme ='test_'
        Mss = 'Mass'
        root = 'Acc_record3/'
        ext = '.txt'
        filen = Nme+str(i)+Mss
        Fn = root+filen+ext
        return Fn      
    if sw==4:
        Nme ='test_'
        Mss = 'Mass'
        root = 'Acc_record4/'
        ext = '.txt'
        filen = Nme+str(i)+Mss
        Fn = root+filen+ext
        return Fn      
    if sw==5:
        Nme ='test_'
        Mss = 'Mass'
        root = 'Acc_record5/'
        ext = '.txt'
        filen = Nme+str(i)+Mss
        Fn = root+filen+ext
        return Fn  
    if sw==6:
        Nme ='test_'
        Mss = 'Mass'
        root = 'Excitation/'
        ext = '.txt'
        filen = Nme+str(i)+Mss
        Fn = root+filen+ext
        return Fn 
#---------------------------- 16.Strain Damage E reduction -------------------#
def Strain_Damage_E_Reductions(Eo,e):
    Esat =  0.817*Eo
    zeta = 46.4
    E =Eo-(Eo-Esat)*(1-np.exp(-zeta*e))
    return E

#---------------------------- 17. Load acceleration data ---------------------#  
# def accelerations(Fn,fi,fj,No,Nf):
#     Acco = np.loadtxt(Fn)
#     Acc = Acco[No:Nf,:]
#     print(len(Acc[:,0]/2))
#     Nc = len(Acc[0,:])
#     Yxx,freq_id,N = BS.PSD_FORMAT(Acc,250,fi,fj)

#     return Yxx,freq_id,Nc,N

#---------------------------- 18. ploting ------------------------------------#  
def ploting_results(freq_id,Yxx):
        plt.figure()
       
        plt.plot(freq_id,Yxx,color = 'blue')
        # plt.plot(freq_id,Yxx)
        plt.yscale('log')
        plt.xlabel('frequency [Hz]')
        plt.ylabel('mm/s^2/Hz^0.5')
        plt.winter()
        plt.show() 
        
#---------------------------- 19. Calling Results ----------------------------#  
def procesing_Data(Cases,i):
    Fn = filecreations(Cases,i,1)
    RNAME = Fn
    R = np.loadtxt(RNAME,delimiter= ' ')
    return R

#---------------------------- 20. Ploting few stuff --------------------------#  
# def visual_stufff(sw,k):
#     Fn = filecreations(sw,k,1)
#     Yxx,freq_id,Nc,N = accelerations(Fn,0.01,50)
#     ploting_results(freq_id,Yxx)
#---------------------------- 21. Mass d4.istribution --------------------------#  
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
        M  = masa[jj,1]+Mb/g
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

        
#---------------------------- 23. Plot Modeshape -----------------------------#      
def plt_modeshape(modid):
    try:
        opsv.plot_mode_shape(modid)
    except:
        ax = plt.gca(projection='3d')
        ax.grid(False)
        
        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')
    
        plt.show()   

#---------------------------- 23. Most unusefull shit   ----------------------#    

def fin_nodes(locs,Const,i):

    if len(locs)>0:
        for k in range(len(locs)):
            print('aca entro con :', k,' con locs ', locs[k])
            bid = np.where(locs[k]== Const)[0][:]
            if len(bid)>0: 
                ll = ops.nodeCoord(Const[i])
                lj = ops.nodeCoord( Const[bid])
                
                
                print(locs[k],ll)
                print(bid,lj)
                plt.plot([ll[0],lj[0]],[ll[1],lj[1]])
                print('este valor esta en el constrain',locs[k])
            else:
                print('iteracion:', k,' este valor no esta en el constrain ', locs[k])
                
    return       
def plt_2D_Modeshapes(Const,conect):
    plt.figure()
    for i in range(len(Const)):
        print('*******************')
        locs =conect[ np.where(Const[i]== conect[:,0])[0][:],1]
        print(locs,len(locs))
        locs1 = conect[np.where(Const[i]== conect[:,1])[0][:],0]
        print(locs1,len(locs1))       
        fin_nodes(locs,Const,i)
        fin_nodes(locs1,Const,i)
        
                            

#---------------------------- 24. plot modeshapes       ----------------------#  
def plot_2d_top_bin():
    plt.figure()
    Const2 = [155, 212, 276, 271, 266, 238, 178, 145, 150,155]
    Const1 = [17,70,138,133,128,97,37,7,12,17]
    for  i in range(len(Const1)-1):
        I = ops.nodeCoord(Const1[i])
        J = ops.nodeCoord(Const1[i+1])
        plt.plot([I[0],J[0]],[I[1],J[1]],Color = 'Blue')
    for  i in range(len(Const2)-1):
        I = ops.nodeCoord(Const2[i])
        J = ops.nodeCoord(Const2[i+1])
        plt.plot([I[0],J[0]],[I[1],J[1]],Color= 'blue')
        plt.ylim([0, 25000])         
#---------------------------- 24. drawing modeshapes   ----------------------#      
def drawing_Modeshapes(U_ops,fs,NM):
    Const2 = [155, 212, 276, 271, 266, 238, 178, 145, 150,155]
    Const1 = [17,70,138,133,128,97,37,7,12,17]
    for i in range(len(Const2)-1):
        I = ops.nodeCoord(Const2[i])
        J = ops.nodeCoord(Const2[i+1])
        plt.plot([I[0],J[0]],[I[1],J[1]],Color = 'Blue')
    for i in range(len(Const1)-1):
        if i < 8:
            modxi = U_ops[i,NM]*fs
            modyi = U_ops[i+18,NM]*fs
            modxj = U_ops[i+1,NM]*fs
            modyj = U_ops[i+19,NM]*fs
        else:
            modxi = U_ops[8,NM]*fs
            modyi = U_ops[8+18,NM]*fs
            modxj = U_ops[0+1,NM]*fs
            modyj = U_ops[0+18,NM]*fs            
        
        I = ops.nodeCoord(Const1[i])
        J = ops.nodeCoord(Const1[i+1])
        plt.plot([I[0]+modxi,J[0]+modxj],[I[1]+modyi,J[1]+modyj],Color= 'blue')
        plt.scatter([I[0]+modxi,J[0]+modxj],[I[1]+modyi,J[1]+modyj],Color= 'red')
        plt.ylim([0, 25000])
    c = 0
    for  i in range(len(Const1)-1,2*len(Const1)-2):
        if i < 2*len(Const1)-1:
            modxi = U_ops[i,NM]*fs
            modyi = U_ops[i+18,NM]*fs
            modxj = U_ops[i+1,NM]*fs
            modyj = U_ops[i+19,NM]*fs
        else:
            modxi = U_ops[2*len(Const1)-2,NM]*fs
            modyi = U_ops[2*len(Const1)-2+18,NM]*fs
            modxj = U_ops[len(Const1)+1,NM]*fs
            modyj = U_ops[len(Const1)+18,NM]*fs            
             
        I = ops.nodeCoord(Const2[c])
        J = ops.nodeCoord(Const2[c+1])
        plt.plot([I[0]+modxi,J[0]+modxj],[I[1]+modyi,J[1]+modyj],Color= 'blue')
        plt.scatter([I[0]+modxi,J[0]+modxj],[I[1]+modyi,J[1]+modyj],Color= 'red')
        plt.ylim([0, 25000])
        c = c+1
        
        
#---------------------------- 25. drawing modeshapes   ----------------------#      
def record_Temperature_Mass(Mb,Mbi,T,E,Mt,Tt,Et,freq,F_M,i):
    freq = np.array(freq)**0.5/(2*np.pi)
    Mt= np.append(Mt,Mb+Mbi*18)
    Tt = np.append(Tt,T)
    Et = np.append(E,Et)
    
    F_M[i,:]= freq
    return Mt,Tt,Et,F_M
    
    
#---------------------------- 26. Save Results   ----------------------#      
 
def writer_txt(sw,Mt,Tt,Et,F_M):
    fn = filecreations(sw,7899871776,0)
    S = np.column_stack((np.round(Mt,4),np.round(Tt,3),np.round(Et,3) ,np.round(F_M,4)))
 
    np.savetxt(fn,S)
    return fn
#---------------------------- 27. Save Results   ----------------------#      

def scatter_3d(fn):
    Data = np.loadtxt(fn)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Data[:,0],Data[:,1],Data[:,3],color = 'blue')
    ax.set_xlabel('Mass [N]')
    ax.set_ylabel('Temperature[C]')
    ax.set_zlabel('Frequency [Hz]')
    
    #---------------------------- 28. F_t_m ---------------------------------#      

def F_t_m(F,M,T,ao,a1,a2,a3,a4,a5):
    f = (ao+a1*M**-a3+a2*T**-a4+a5*M)-F
    return f

def F_t_m3(F,M,T,ao,a1,a2,a3,a4,a5):
    f = ao+a1*M**-a3+a2*T**-a4+a5*M
    return f

def fit_model_t_mass(fn,mod):
    from scipy.optimize import least_squares
    Data = np.loadtxt(fn)
    M = Data[:,0]
    T = Data[:,1]
    F = Data[:,2+mod]

    fun = lambda x: F_t_m(F,M,T,x[0],x[1],x[2],x[3],x[4],x[5])
    x0 = [4,1,1,0.5,0.6,1] 
    opt = least_squares(fun, x0, loss='soft_l1', f_scale=0.1)
    return opt
    
    #---------------------------- 29. scatter -------------------------------#      
    
def scatter_3d_fit(fn,xpt,mod):
    Data = np.loadtxt(fn)
    fig = plt.figure()
    Z = F_t_m3(Data[:,2+mod],Data[:,0],Data[:,1],xpt[0],xpt[1],xpt[2],xpt[3],xpt[4],xpt[5])
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Data[:,0],Data[:,1],Data[:,2+mod],color = 'blue')
    ax.scatter(Data[:,0],Data[:,1],Z,color = 'red')
    ax.set_xlabel('Mass [N]')
    ax.set_ylabel('Temperature [C]')
    ax.set_zlabel('Frequency [Hz]')

    #---------------------------- 30. Surface PLot ---------------------------#      

def surface_model(fn,xpt,mod):
    from matplotlib.colors import LightSource
    plt.close('all')
    from matplotlib import cm
    Data = np.loadtxt(fn)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    X = np.arange( 108333, 3.5*10**6, 10000)
    Y = np.arange(20 , 60, (60-20)/len(X))
    X, Y = np.meshgrid(X, Y)
    Z = F_t_m3(Data[:,2+mod],X,Y,xpt[0],xpt[1],xpt[2],xpt[3],xpt[4],xpt[5])
    ls = LightSource(270, 45)
    rgb = ls.shade(Z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
    surf = ax.plot_surface(X, Y, Z,  rstride=1, cstride=1, facecolors=rgb,
                   linewidth=0, antialiased=False, shade=False)
    ax.scatter(Data[:,0],Data[:,1],Data[:,2+mod],s = 1,color = 'green', marker='x')
    ax.set_xlabel("Mass [N]"); ax.set_ylabel("Temperature [C]"); ax.set_zlabel("Frequency [Hz]")
     
    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)

     
#---------------------------- 28. Uncertainties  ---------------------------------#    
  
def unceratinties_model(fn,xpt,mod):
    Data = np.loadtxt(fn)
    M = Data[:,0]
    T = Data[:,1]
    F = Data[:,2+mod]
    from scipy import optimize
    # cost function
    def like_regretion(F,M,T,xo):
        er = np.sum(-np.log(stats.norm.pdf( F_t_m3(F,M,T,xo[0],xo[1],xo[2],xo[3],xo[4],xo[5]),F,xo[6])))
        return er 
    Err = lambda x:  like_regretion(F,M,T,x)
    xo = np.append(xpt,0.01)
    xptf = optimize.fmin(Err,xo)
    h = np.random.normal( F_t_m3(1,2*10**6,30,xo[0],xo[1],xo[2],xo[3],xo[4],xo[5]),xptf[-1]/2,1000)
    plt.figure()
    plt.hist(h,bins = 20)
    plt.title('Bivariate Regretation')
    return xptf



#---------------------------- 29. Crusher Excitation -------------------------#    
    
    
def Timehistory_crusher(Fn,Ndata,Acc_Series):
    
    Nrsp = [155, 212, 276, 271, 266, 238, 178, 145, 150,  17,  70, 138, 133,\
    128,  97,  37,   7,  12]
    ops.wipeAnalysis()
    freq = ops.eigen(4)
    ops.timeSeries('Path', 1, '-dt', 1/500, '-filePath',Acc_Series, '-factor',1)
    ops.recorder('Node', '-file', Fn  ,'-timeSeries', 1,'-dT', 1/250,'-node', *Nrsp,'-dof',1,'accel')
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)    
    dampRatio = 0.002
    ops.rayleigh(0., 0., 0., 2*dampRatio/freq[0])
    # create the analysis
    ops.constraints(' Transformation')
    ops.numberer('Plain')    # renumber dof's to minimize band-width (optimization), if you want to
    ops.system('BandGeneral') # how to store and solve the system of equations in the analysis
    ops.algorithm('Linear')     # use Linear algorithm for linear analysis
    ops.integrator('Newmark', 0.5, 0.25)    # determine the next time step for an analysis
    ops.analysis('Transient')   # define type of analysis: time-dependent
    ops.analyze(int(Ndata/2), 1/250)     # apply time steps in analysis
    ops.remove('recorders')    
#------------------------------- 32. Conenction curve -------------------------------#    

def curve_stiffness(Eo,Fy,eEo,ef,M):

    R = 10
    E1 = (Fy-(Fy*2/3))/(ef-eEo)
    b = E1/Eo
    e = np.linspace(0,ef,100)
    sigma = b*e + (1-b)*e/(1+e**R)**(1/R)
    plt.plot(e,sigma*Fy,color = 'blue',)
    if M =='M':
        plt.xlabel('mrad [$\\theta$] ')
        plt.ylabel('M [kMm]')
        plt.legend(['Menegotto-Pinto'])

#------------------------------- 33. Conenction curve ------------------------#    
    
    
def ploting3d(Nodes,conect,sw,ax,c,alfa):
    # from matplotlib.widgets import Button
    if ax == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.subplots_adjust(bottom=0.2)
    

    for i in range(len(conect[:,0])):
        ii = conect[i,0]
        jj = conect[i,1]
        corrii = Nodes[int(ii-1),:]
        corrjj = Nodes[int(jj-1),:]
        ax.plot([corrii[0],corrjj[0]],[corrii[1],corrjj[1]],[corrii[2],corrjj[2]],color = c)
        
        
                   
    # Get rid of colored axes planes
    # First remove fill
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    # Now set color to white (or whatever is "invisible")
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    
    # Bonus: To get rid of the grid as well:
    ax.grid(False)
    
    ax.scatter(Nodes[:,0],Nodes[:,1],Nodes[:,2],color= 'red',s= 2)
    ax.set_zlim([0,18000])
    ax.set_ylim([0,25000])
    ax.set_xlim([0,25000])
    ax.view_init(90,-90 )
    ax.set_axis_off()
      
    
#------------------------------- 34. Conenction curve -----------------------#    
def subploting(Nodes,conect):
    fig = plt.figure()
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    ploting3d(Nodes,conect,1,ax)
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    ploting3d(Nodes,conect,1,ax)  
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    ploting3d(Nodes,conect,1,ax)       
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    ploting3d(Nodes,conect,1,ax)
    
    
#------------------------------- 35. Get Modeshapes --------------------------#    

def getmodeshapes(Nodes,Nm):
    Nd = Nodes.shape[0]
    Modes = np.zeros([Nd,3])
    for i in range(Nd):
        Modes[i,0]=ops.nodeEigenvector(int(i+1), Nm)[0]
        Modes[i,1]=ops.nodeEigenvector(int(i+1), Nm)[1]
        Modes[i,2]=ops.nodeEigenvector(int(i+1), Nm)[2]
        
    return Modes

#------------------------------- 36. subplotng modeshape ---------------------#    
def subploting_mod(Nodes,conect,FS):
    breakpoint()
    fig = plt.figure()
    #-------- Mod1 ------------------------------#
    ax = fig.add_subplot(2, 3, 1, projection='3d')
    Modes = getmodeshapes(Nodes,1)
    Nodes1 = Nodes + FS*Modes
    ploting3d(Nodes1,conect,1,ax,'blue',1)
    ploting3d(Nodes,conect,1,ax,'gray',0.1)
    #-------- Mod2 ------------------------------#
    ax = fig.add_subplot(2, 3, 2, projection='3d')
    Modes = getmodeshapes(Nodes,2)
    Nodes2 = Nodes + FS*Modes
    ploting3d(Nodes2,conect,1,ax,'blue',1)
    ploting3d(Nodes,conect,1,ax,'gray',0.1)
    #-------- Mod3 ------------------------------#
    ax = fig.add_subplot(2, 3, 3, projection='3d')
    Modes = getmodeshapes(Nodes,3)
    Nodes3 = Nodes + FS*Modes
    ploting3d(Nodes3,conect,1,ax,'blue',1)
    ploting3d(Nodes,conect,1,ax,'gray',0.1)
    #-------- Mod4 ------------------------------#
    ax = fig.add_subplot(2, 3, 4, projection='3d')
    Modes = getmodeshapes(Nodes,4)
    Nodes4 = Nodes + FS*Modes
    ploting3d(Nodes4,conect,1,ax,'blue',1)
    ploting3d(Nodes,conect,1,ax,'gray',0.1)
    #-------- Mod5 ------------------------------#
    ax = fig.add_subplot(2, 3, 5, projection='3d')
    Modes = getmodeshapes(Nodes,5)
    Nodes4 = Nodes + FS*Modes
    ploting3d(Nodes4,conect,1,ax,'blue',1)
    ploting3d(Nodes,conect,1,ax,'gray',0.1)
    #-------- Mod6 ------------------------------#
    ax = fig.add_subplot(2, 3, 6, projection='3d')
    Modes = getmodeshapes(Nodes,6)
    Nodes4 = Nodes + FS*Modes
    ploting3d(Nodes4,conect,1,ax,'blue',1)
    ploting3d(Nodes,conect,1,ax,'gray',0.1)
    
        
def desing_ratio(Nodes,conect,ax):
    import matplotlib.colors
    cmap = plt.cm.rainbow
    code = np.random.uniform(0,1,509)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
    
    
    if ax == 0:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.subplots_adjust(bottom=0.2)
    

    for i in range(len(conect[:,0])):
        ii = conect[i,0]
        jj = conect[i,1]
        corrii = Nodes[int(ii-1),:]
        corrjj = Nodes[int(jj-1),:]
        ax.plot([corrii[0],corrjj[0]],[corrii[1],corrjj[1]],[corrii[2],corrjj[2]],color=cmap(norm(code[i])))
        if code[i] > 0.7:
            xtx = corrii[0]+(-corrii[0]+corrjj[0])/2
            ytx = corrii[1]+(-corrii[1]+corrjj[1])/2
            ztx = corrii[2]+(-corrii[2]+corrjj[2])/2
            ax.text(xtx, ytx,ztx,str(round(code[i],1)), style='italic',fontsize=5)
    #-------- Cleaning axis -----------------------#
        
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    
    ax.grid(False)
    #-------- Scatter -----------------------#
    ax.scatter(Nodes[:,0],Nodes[:,1],Nodes[:,2],color= 'red',s= 2)
    #-------- Limit axis --------------------#
    ax.set_zlim([0,18000])
    ax.set_ylim([0,25000])
    ax.set_xlim([0,25000])
    ax.view_init(90,-90 )
    ax.set_axis_off() 
    
    #-------- side bar :) -----------------------#
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([]) 
    fig.colorbar(sm)
    plt.title('Design Ratio')
    plt.show()