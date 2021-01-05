#                               SSI COV 
#                               Functions

import numpy as np
import matplotlib.pyplot as plt 
from scipy.cluster.hierarchy import linkage,fcluster,fclusterdata
plt.close('all')
# --------------------------- 1. Load Data -----------------------------------#
def load_data(fn,delim):
    try:
        Acc= np.loadtxt(fn,delimiter = delim)

        Nc  = len(Acc[0,:]) 
        N   = len(Acc[:,0])
        return Acc,Nc,N
    except:
        print('error: could not load the txt file')        
        return [],[],[]
    
# --------------------------- 2. pot_Acc -------------------------------------#

def plot_ACC(Acc,io,ie):
    plt.figure()
    for i in range(io,ie):
        plt.plot(Acc[:,i])
        plt.xlabel('Numner of data Points')
        plt.ylabel('Acceleration ')
        plt.title('Acceleration Record')
        plt.show()


# --------------------------- 3. NexT ----------------------------------------#
def NexT(acc,fs,Ts,Nc):
    dt = 1/fs
    M = round(Ts/dt)

    IRF = np.zeros((Nc,Nc,M-1),dtype = complex) 
    for oo in range(Nc):
        for jj in range(Nc):
            y1 = np.fft.fft(acc[:,oo])
            y2 = np.fft.fft(acc[:,jj])
            h0 = np.fft.ifft(y1*y2.conj())
            #plt.plot(np.real(h0[0:M]))
            IRF[oo,jj,:] = np.real(h0[0:M-1])
           

    t = np.linspace(0,1,M-1)*dt
    if Nc ==1:
        IRF = np.squeeze(IRF)
        IRF = IRF/IRF[0]
    
    return IRF
    
    
# --------------------------- 4. blockToeplitz -------------------------------#
def blockToeplitz(h):
    
    N1 = round(h.shape[2]/2)-1
    M = h.shape[1]
    
    # print('M:     ',M)
    # print('Size T1',(N1+1)*M-1, 's ' , (N1+1)*M-1)
    T1= np.zeros(((N1)*M,(N1)*M))



    for oo in range(N1):
        for ll in range(N1):
            
            # print((oo+1)*M-1,(ll+1)*M-1,len(h[:,:,N1+oo-ll]))
            T1[(oo)*M:(oo+1)*M,(ll)*M:(ll+1)*M] = h[:,:,N1-1+oo-ll+1]
    
    [U,S,V] = np.linalg.svd(T1)
    # print(U[:,0])
    return U,S,V,T1
# --------------------------- 5. Modal id ------------------------------------#

def modalID(U,S,Nmodes,Nyy,fs):
    S = np.diag(S)
    dt = 1/fs
    O = np.matmul(U[:,0:Nmodes],np.sqrt(S[0:Nmodes,0:Nmodes]))

    IndO = min(Nyy,len(O[:,0]))

    C = O[0:IndO,:]

    jb =O.shape[0]/IndO
    ao = int((IndO)*(jb-1))
    bo = int(len(O[:,0])-(IndO)*(jb-1))
    co = len(O[:,0])
    A =np.matmul( np.linalg.pinv(O[0:ao,:]),O[bo:co,:])
    [Vi,Di] = np.linalg.eig(A)
    mu = np.log(np.diag(np.diag(Vi)))/dt
    fno = np.abs(mu)/(2*np.pi)
    fn = fno[np.ix_(*[range(0,i,2) for i in fno.shape])]
    zetaoo = -np.real(mu)/np.abs(mu)
    zeta =  zetaoo[np.ix_(*[range(0,i,2) for i in zetaoo.shape])]
    phi0 = np.real(np.matmul(C[0:IndO,:],Di))

    # phi = phi0[np.ix_(*[range(0,i,2) for i in phi0.shape])]
    phi = phi0[:,1:-1:1]


    return fn,zeta,phi

    """
    multiplicacion de matrices
    numeros imaginarios
    cosa que comienzas en 2 y en 1 
    """
# --------------------------- 6. Stability Check  ----------------------------#

def  stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1):
    eps_freq = 0.1 
    eps_zeta = 0.4 
    eps_MAC = 0.05
    stablity_status = np.array([])
    fn = np.array([])
    zeta = np.array([])
    phi = []
    MAC= np.array([])
    stabStatus = 0
     # frequency stability
    N0 = len(fn0);
    N1 = len(fn1);
    for i in range(N0):
        for j in range(N1):
            stab_fn = errorcheck(fn0[i],fn1[j],eps_freq)
            stab_zeta = errorcheck(zeta0[i],zeta1[j],eps_zeta)
            try:
                stab_phi,dummyMAC = getMAC(phi0[:,i],phi1[:,j],eps_MAC)
           
            except:
                print(i)
                print(j)
                
            if stab_fn==0:
                stabStatus =0
            elif stab_fn == 1 & stab_phi == 1 & stab_zeta == 1:
                stabStatus = 1
            elif stab_fn == 1 & stab_zeta ==0 & stab_phi == 1:
                stabStatus = 2
            elif stab_fn == 1 & stab_zeta == 1 & stab_phi ==0:
                stabStatus = 3
            elif stab_fn == 1 & stab_zeta ==0 & stab_phi ==0:
                stabStatus = 4
    
            
         
            fn = np.append(fn,fn1[j])
          
            zeta = np.append(zeta,zeta1[j])
            dummyy = phi1[:,j]
            phi.append(dummyy)
          
            MAC= np.append(MAC,dummyMAC)
            try:
                
                stablity_status = np.append(stablity_status,stabStatus)
            except:
                # print('stab_fn : ',stab_fn)
                # print('stab_phi ',stab_phi)
                breakpoint()
                if stab_fn==0:
                    stabStatus =0
                elif stab_fn == 1 & stab_phi == 1 & stab_zeta == 1:
                    stabStatus = 1
                elif stab_fn == 1 & stab_zeta ==0 & stab_phi == 1:
                    stabStatus = 2
                elif stab_fn == 1 & stab_zeta == 1 & stab_phi ==0:
                    stabStatus = 3
                elif stab_fn == 1 & stab_zeta ==0 & stab_phi ==0:
                    stabStatus = 4
                    
                stablity_status = np.append(stablity_status,stabStatus)

    phi = np.asarray(phi )
    phi = phi.T

    iddd = np.argsort(fn,-1)
   
    fn = fn[iddd]
    zeta = zeta[iddd]
    phi = phi[:,iddd]
    MAC = MAC[iddd];
    stablity_status = stablity_status[iddd];
        
    
    return fn,zeta,phi,MAC,stablity_status

# --------------------------- 7. error check ---------------------------------#
def errorcheck(xo,x1,eps):
    if abs(1-xo/x1)<eps:
        y = 1
    else:
        y = 0
    return y

# --------------------------- 8. Get Mac -------------------------------------#

def getMAC(x0,x1,eps):
    Num = np.abs(np.matmul(x0[:],(x1[:].reshape(-1,1))))**2
     
    # print('¨¨¨¨¨¨¨¨')
    # print(Num)
    # print('¨¨¨¨¨¨¨¨')
    D1 = np.matmul(x0[:],(x0[:].reshape(-1,1)))
    D2 = np.matmul(x1[:],(x1[:].reshape(-1,1)))
    dummpyMac = Num/(D1*D2)
    # print(dummpyMac)
    if dummpyMac >(1-eps):
        y = 1
    else:
        y = 0 
    return  y,dummpyMac
# --------------------------- 9. flip dictionary -----------------------------#

def flip_dic(a):
    
    from collections import OrderedDict
    d = OrderedDict(a)
    dreversed = OrderedDict()
    for k in reversed(d):
        dreversed[k] = d[k]
        
    return dreversed


# --------------------------- 9.Get stable Poles -----------------------------#

def getStablePoles(fn,zeta,phi,MAC,stablity_status):
    fnS = np.array([])
    zetaS = np.array([])
    phiS = []
    MACS= np.array([])

    for i in range(len(fn)):
        for j in range(len(stablity_status[i])):
            if stablity_status[i][j]==1:

                fnS = np.append(fnS,fn[i][j])
                zetaS = np.append(zetaS,zeta[i][j])
                dummyyS= phi[i][:,j]
                phiS.append(dummyyS)
                MACS= np.append(MACS,MAC[i][j])
                
    phiS = np.asarray(phiS )
 
    phiS = phiS.T
    #  remove negative damping
    fnS = np.delete(fnS,np.where(zetaS<0))
    phiS= np.delete(phiS,np.where(zetaS<0),1)
    MACS =  np.delete(MACS,np.where(zetaS<0))
    zetaS = np.delete(zetaS,np.where(zetaS<0))
    
    for oo in range(phiS.shape[1]):
        phiS[:,oo] = phiS[:,oo]/np.linalg.norm(phiS[:,oo])
        if np.diff(phiS[0:1,oo]) < 0:
            phiS[:,0] = -phiS[:,oo]
            
    return fnS,zetaS,phiS,MACS

# --------------------------- 10.Plot stabilization  --------------------------#

def plot_stabilization_diagram_2(data):

    plt.figure(figsize=(8, 6), dpi=80)
    plt.subplot(1, 1, 1)

    for pos, row in enumerate(data):
        for pos2, item in enumerate(row):
            plt.plot(pos, item, 'o')

    plt.show()
 
# --------------------------- 11.Cluster  --------------------------#
def ClusterFun(fn0,zeta0,phi0): 

    Nsamples = phi0.shape[1]
    # print(Nsamples)
    pos = np.zeros((Nsamples,Nsamples))
    for i in range(Nsamples):
        for j in range(Nsamples):
             stab_phi,MAC0 = getMAC(phi0[:,i],phi0[:,j],0.05)
             pos[i,j] = np.abs((fn0[i]-fn0[j])/(fn0[j]))
    
    Z = linkage(pos,'single','euclidean')
    myClus = fcluster(Z,0.1,criterion = 'distance')

    Ncluster = max(myClus)
    ss= 0
    fn = {}
    zeta = {}
    phi = {}

    for rr in range(Ncluster):
        # print('Iteracion-  ',rr , 'len')
        if len(myClus[np.where(myClus == rr)])>5:
         
            dummyZeta = zeta0[np.where(myClus==rr)]
            dummyFn = fn0[np.where(myClus==rr)]
            dummyPhi = phi0[:,np.where(myClus==rr)[0]]
            valMin = max(0,(np.quantile(dummyZeta,0.25)-abs(np.quantile(dummyZeta,0.75)-np.quantile(dummyZeta,0.25))*1.5))
            valMax = np.quantile(dummyZeta,0.75)+abs(np.quantile(dummyZeta,0.75)-np.quantile(dummyZeta,0.25))*1.5
            
            dummyFn = np.delete(dummyFn,np.where((dummyZeta>valMax)|(dummyZeta<valMin)))
            dummyPhi = np.delete(dummyPhi ,np.where((dummyZeta>valMax)|(dummyZeta<valMin))[0],1)
            dummyZeta = np.delete(dummyZeta,np.where((dummyZeta>valMax)|(dummyZeta<valMin)))
            
      
            fn[ss]= dummyFn
            zeta[ss] = dummyZeta
            
            phi[ss] = dummyPhi
            ss= ss +1 
   
    return fn,zeta,phi

# --------------------------- 12. Results  --------------------------#
def Cluster_Resuls(fn,zeta):
    std,mean,mean_d,std_d = np.array([]),np.array([]),np.array([]),np.array([])
    for i in range(len(fn)):
        std = np.append(std,np.std(fn[i])/len(fn[i]))
        mean = np.append(mean,np.mean(fn[i]))
        mean_d = np.append(mean_d,np.mean(zeta[i]))
        std_d = np.append(std_d,np.std(zeta[i]))
    
    
    
    
    idd = np.argsort(std)
    std = std[idd]
    mean = mean[idd]
    mean_d = mean_d[idd]
    std_d = std_d[idd]
    
    
    # print('Identified parameters: ')
    # print('Modal Frequencies: ')
    # print(mean)
    # print('Standard Deviation ')
    # print(std)
    # print('Damping')
    # print(mean_d)
    # print('Standard Deviation ')
    # print(std_d)
    
    return mean,mean_d
                
# --------------------------- 13. Results  --------------------------#

def shape_modeshape(phi,fn):
    breakpoint()
    clases = len(phi)
    phi_shape = {}
    for i in range(clases):
        
        col= len(fn[i])
        row =int(len(phi[i])/col)
        phi_shape[i] = phi[i].reshape((row,col))

                
            
            
        
        
            

        




































    
