from pathlib import Path

import numpy as np
import openseespy.opensees as ops


current_path = Path(str(__file__)).parent

#---------------------------- 1. Topology  -----------------------------------#
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
    Names = np.genfromtxt(current_path/'Nmes.txt',dtype='str',delimiter = '')
    return Nodes, conect, idele,Members,Names
#---------------------------- 2. Hailcreek -----------------------------------#
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
#---------------------------- 3. Internal Forces -----------------------------#
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
#---------------------------- 5. Mass distribution ---------------------------#  
def mass_distirbution(masa,Mbo,MB1,MB2,dist_b1,dist_b2,LfM):

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
        ops.load(masa[jj,0],0,0,LfM*M*g,0,0,0)
    for i in range(len(Const2)):
        Mb = Mbo+MB2*dist_b2[i]
        jj = int(Const2[i])-1
        M  = masa[jj,1]+Mb/g
        ops.mass(masa[jj,0],M,M,M,0,0,0)
        ops.load(masa[jj,0],0,0,LfM*M*g,0,0,0)
#---------------------------- 6. Mass Shape Distribution ---------------------#     
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
#---------------------------- 7. zero length ---------------------------------#
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
#---------------------------- 8. Load Cases ----------------------------------#
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
        loads = np.loadtxt(current_path/'Loads_wx.txt')
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
                    # print(tag)
                    tag = loads[i,0]
                    ter = np.sum(xo)/np.sum(zi)
                    wi = loads[i,1]*ter
                    # print('¨*_****',wi,'****')
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
        loads = np.loadtxt(current_path/'Loads_wy.txt')
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
#---------------------------- 9. Desing Stress Forces ------------------------#
def Desing_Stress_Forces(ld,conect,idele,Ew):
    Ful_Sec = np.loadtxt(current_path/'Members_full.txt')
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
                np.savetxt(current_path/'2.txt',s,delimiter = '  ')
    else:       
          if ld >= 5:
              return 
          else:
              if ld == 3:
                  Loads = np.loadtxt(current_path/'Loads_wx.txt')   
              else:
                 Loads = np.loadtxt(current_path/'Loads_wy.txt') 
                      
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
                    np.savetxt(current_path/'496.txt',s,delimiter = '  ')               
#---------------------------- 10. Forces -------------------------------------#
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
#---------------------------- 11.Rot Tranf -----------------------------------#
def rot_transf_3d(ex, ey, ez, g):

    Lxyz = np.array([ex[1]-ex[0], ey[1]-ey[0], ez[1]-ez[0]])

    L = np.sqrt(Lxyz @ Lxyz)

    z = np.zeros((3, 3))

    G = np.block([[g, z, z, z], 
                  [z, g, z, z],
                  [z, z, g, z],
                  [z, z, z, g]])

    return G, L
#---------------------------- 12. Internal Forces-----------------------------#
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

def all_loads(idpattern,conect,Members,gama,idele,LF,LfM):
    """
    first iterate over all the elements
    """
    ops.wipeAnalysis()
    ops.timeSeries('Linear',2)
    ops.pattern("Plain", idpattern, 2)
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
        m =Members[int(idele[i])-1,3]*gama*LfM
        
        if np.linalg.norm(vecxz) == 0:
            # ops.eleLoad('-ele',i+1,'-type','-beamUniform',0, m,0) 
            
            Ew[i+1] = {i+1: ['-beamUniform',0, 0,m]}
            
            
        else:
            # ops.eleLoad('-ele',i+1,'  -type','-beamUniform',m,0,0)  
            Ew[i+1] = {i+1: ['-beamUniform',m, 0,0]}
        # breakpoint()
    Ew = windloadcases(Ew,LF)
    # breakpoint()
    return Ew
        
            
            
def windloadcases(Ew,LF):
        
        loads = np.loadtxt(current_path/'Loads_wx.txt')

        xo = [1,0,0]
        Ldfactor = np.matmul(LF,xo)
        Ew = winds_loads(Ew,loads,xo,Ldfactor)
        #Pointloads wx
        Ploads = np.loadtxt(current_path/'Pointloadwx_b1.txt',delimiter = ',')
        assign_poinloads(Ldfactor,Ploads)
        Ploads = np.loadtxt(current_path/'Pointloadwx_b2.txt',delimiter = ',')
        assign_poinloads(Ldfactor,Ploads)


        
        loads = np.loadtxt(current_path/'Loads_wy.txt')
        xo = [0,1,0]
        Ldfactor = np.matmul(LF,xo)
        Ew = winds_loads(Ew,loads,xo,Ldfactor)
        #Pointloads wy
        Ploads = np.loadtxt(current_path/'Pointloadwy_b1.txt',delimiter = ',')
        assign_poinloads(Ldfactor,Ploads)
        Ploads = np.loadtxt(current_path/'Pointloadwy_b2.txt',delimiter = ',')
        assign_poinloads(Ldfactor,Ploads)
        
        return Ew       
        
        # Y = [1,0,0]
        # Z = [0,1,0]
        
def winds_loads(Ew,loads,xo,Lfac):
        # breakpoint()
        z = [0,0,1]
        for i in range(len(loads[:,0])):
            tag = loads[i,0]
            yi = ops.eleResponse(tag, 'ylocal')
            zi = ops.eleResponse(tag, 'zlocal')
            if np.linalg.norm(np.cross(xo,yi)) == 0:
                    ter = np.sum(xo)/np.sum(yi)
                    wi = loads[i,1]*ter
                    # ops.eleLoad('-ele',tag,'-type','-beamUniform', wi,0,0)
                    Ew[int(tag)][int(tag)][2] =+ wi*Lfac
                    # Ew[tag] = {tag: ['-beamUniform',wi,0,0]}
            else:
                    # breakpoint()
                    # print(tag)
             
                    ter = np.sum(xo)/np.sum(zi)
                    wi = loads[i,1]*ter
                    # print('¨*_****',wi,'****')
                    # ops.eleLoad('-ele',tag,'-type','-beamUniform', 0,wi,0)
                    # Ew[tag] = {tag: ['-beamUniform',0, wi,0]}
                    Ew[int(tag)][int(tag)][2] =+ wi*Lfac
        return Ew    
    
def assignloads_dist(Ew):
    for i in range(len(Ew)):
        wi,wj,wk = Ew[int(i+1)][int(i+1)][1:]
        ops.eleLoad('-ele',int(i+1),'-type','-beamUniform', wi,wj,wk)

def assign_poinloads(Ldfactor,Ploads):
    """
    this part i assing the point loads for the wind
    id Vx Vy Vz
    bae
    """
    for i in range(Ploads.shape[0]):
        tag = int(Ploads[i,0])
        Px =  Ploads[i,1]*Ldfactor
        Py =  Ploads[i,2]*Ldfactor
        Pz =  Ploads[i,3]*Ldfactor
        LOADS = [Px,Py,Pz,0,0,0]
        ops.load(tag,*LOADS)

def assign_point_seismic(Lfactor):
    Sl = np.loadtxt(current_path/'Seismic_Load_Hk.txt')
    for i in range(Sl.shape[0]-1):
        Px = Sl[i,0]*Lfactor[0]*1000
        Py = Sl[i,1]*Lfactor[1]*1000
        Pz = Sl[i,2]*Lfactor[2]*1000
        LOADS = [Px,Py,Pz,0,0,0]
        tag = int(i+1)
        ops.load(tag,*LOADS)
        
def assign_point_truck(Lfactor):
    Sl = np.loadtxt(current_path/'Truckload.txt')
    
    for i in range(Sl.shape[0]-1):
        Px = Sl[i,0]*Lfactor
        Py = Sl[i,1]*Lfactor
        Pz = Sl[i,2]*Lfactor
        LOADS = [Px,Py,Pz,0,0,0]
        tag = int(i+1)
        ops.load(tag,*LOADS)    

                            

#---------------------------- 9. Desing Stress Forces ------------------------#
def Desing_Stress_Forces2(ld,conect,idele,Ew,Names):
    

    Svonmisses = {}
    Design_ratio = {}
    Ful_Sec = np.loadtxt(current_path/'Members_full.txt')  
    Cpcity = np.loadtxt(current_path/'capacity.txt',delimiter =',')
            
    for i in range(len(conect[:,1])):
        #(N, Vy, Vz, T, My, Mz)
        
        tag = int(i+1)
        s,Mz =Forces(tag,Ew[tag])
        My = s[:,4]
        Vz = s[:,2]
        Vy = s[:,1]
        N  = s[:,0]
        
        Zy = Ful_Sec[int(idele[int(tag-1)]-1),15]
        Zx = Ful_Sec[int(idele[int(tag-1)]-1),16]
        AA = Ful_Sec[int(idele[int(tag-1)]-1),6]
        Asy = Ful_Sec[int(idele[int(tag-1)]-1),12]
        Asx = Ful_Sec[int(idele[int(tag-1)]-1),11]
        
        Mz_p =max(Mz)/1000000
        Mz_n =min(Mz)/1000000
        My_p =max(My)/1000000
        My_n =min(My)/1000000
        
        Vz_p =max(Vz)
        Vz_n =min(Vz)
        Vy_p =max(Vy)
        Vy_p =min(Vy)
        
        N_p  = max(N)
        N_n  = min(N)
        
        sigMz_p = Mz_p/Zy
        sigMz_n = Mz_n/Zy
        
        sigMy_p = My_p/Zx
        sigMy_n = My_n/Zx
        
        SgN_p = N_p/AA
        SgN_n = N_n/AA
        
        tauz_p = Vz_p/Asy
        tauz_n = Vz_n/Asy
          
        tauy_p = Vy_p/Asx
        tauy_n = Vy_p/Asx
        # breakpoint()
        if  Cpcity[int(idele[i]-1),2] <1:
            dratio  = max(abs(N))/ (Cpcity[int(idele[i]-1),1]*1000)
        else:
            dt_N = max(abs(N))/ (Cpcity[int(idele[i]-1),1]*1000)
            dt_Mz = max(abs(Mz)*0.0000010)/ Cpcity[int(idele[i]-1),2]
            dt_My = max(abs(My)*0.0000010)/ Cpcity[int(idele[i]-1),3]
            Alldratio = [dt_N,dt_Mz,dt_My]
            dratio = max(np.abs(Alldratio))
        
        sig_N =sigMz_p+sigMz_n+sigMy_p+sigMy_n + SgN_p+SgN_n
        tau_T = tauz_p+tauz_n+tauy_p+ tauy_n
        sig_Vommises = (sig_N**2+3*tau_T**2)**0.5
        # print('Element id ',tag,' ',Names[int(idele[i]-1)],' Vonnmisses: ',round(sig_Vommises,3),'Desing ratio',dratio)
        # breakpoint()
        Svonmisses[tag] = [sig_Vommises, dratio,dt_N,dt_Mz,dt_My]
        Design_ratio[tag] = dratio
    return Svonmisses,Design_ratio
  
def run_model():
    ops.constraints('Transformation')
    ops.algorithm('Linear')
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.integrator('LoadControl', 1)
    ops.analysis('Static')
    ops.analyze(1)
    
def Realiability(Dr):
    # breakpoint()
    pf = np.zeros((Dr.shape[0],1))
    for i in range(Dr.shape[0]):
        beta = 5.86*Dr[i]**2-15.83*Dr[i]+9.97
        # breakpoint()
        pf[i,0] = prob(beta)
    return pf
        
def prob(beta):
    
    if 0 < beta < 1.3:	
        prob = -0.690*beta+1.000
    if 1.3 <= beta < 2.3:	
        prob = -0.09*beta+0.217
    if 2.3 <= beta < 3.1:	
        prob = -0.010*beta+0.035
    if 3.1 <= beta < 3.7:	
        prob = -0.0015*beta+0.0057
    if 3.7 <= beta < 4.2:	
        prob = 10**-5
    if 4.2 <= beta < 4.7:	
        prob = 10**-6
    if 4.7 <= beta:	
        beta = 4.7
        prob = 10**-7

    return prob

def Likelyhoood(prob,cases):
    if cases ==1:
        P = 1/50
        Lk = prob*P
        return Lk
    if cases ==2:
        P = 1
        Lk = prob*P
        return Lk
    if cases ==3:
        P = 1/500
        Lk = prob*P
        return Lk
    if cases ==4:
        P = 1/5
        Lk = prob*P
        return Lk
        