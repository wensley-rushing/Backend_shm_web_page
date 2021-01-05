import backend.SSIPascoBridge.SSI.SSI_Backend as SSI 
import matplotlib.pyplot as plt 
import numpy as np 
from scipy import signal


def SSI_COV_AD(Acc,fs,Ts,Nc,Nmax,Nmin):
# --------------------------- 3. NexT ----------------------------------------#
    #breakpoint()
    IRF = SSI.NexT(Acc,fs,Ts,Nc)

    [U,S,V,T] = SSI.blockToeplitz(IRF)

    # plt.figure()

        
    # fig, ax1 = plt.subplots()
    # color = 'blue'
    # for i in range(Nc):
    #     # freqss,pxx = signal.welch(Acc[:,i],fs,nfft = int(len(Acc[:,i])/2)+1,nperseg=int(len(Acc[:,i])/2))
    #     ax1.set_xlabel('frequency [Hz]')
    #     ax1.set_ylabel('Power Spectral Density')
    
    #     # ax1.plot(freqss, 10*np.log(pxx))
    #     plt.xlim([0,int(fs/2)+1])
            
    #     ax1.tick_params(axis='y', labelcolor=color)
    # ax2 = ax1.twinx() 
    # color = 'red'
    # ax2.set_ylabel('order', color=color) 
    
    
    
    kk=0
    
    # print('******* Identified poles **********')
    fn2,zeta2,phi2,MAC,stablity_status = {},{},{},{},{}
    for i in range(Nmax,Nmin,-1):
        if kk == 0:
            fn0,zeta0,phi0 = SSI.modalID(U,S,i,Nc,fs)
    
        else:
            fn1,zeta1,phi1 = SSI.modalID(U,S,i,Nc,fs)
            
            # print(fn1)
            # ax2.plot(fn1,i*np.ones(len(fn1)),'x',color=color)
            # ax2.tick_params(axis='y', labelcolor=color)        
            [a,b,c,d,e] = SSI.stabilityCheck(fn0,zeta0,phi0,fn1,zeta1,phi1)
        
      
            fn2[kk-1]=a
            zeta2[kk-1]=b
            phi2[kk-1]=c
            MAC[kk-1]=d
            stablity_status[kk-1]=e
            
            fn0=fn1
            zeta0=zeta1
            phi0=phi1  
            
        kk = kk +1
    
    # --------------------------- 9. flip dictionary -----------------------------#
    fn2 = SSI.flip_dic(fn2)
    zeta2 = SSI.flip_dic(zeta2)
    phi2 = SSI.flip_dic(phi2)
  
    
    fnS,zetaS,phiS,MACS = SSI.getStablePoles(fn2,zeta2,phi2,MAC,stablity_status)

    
    fn,zeta,phi = SSI.ClusterFun(fnS,zetaS,phiS)
    
    # print('********************')
    # print('******** cluster :D *********')
    # print(fn)
    # print('********  cluster :D  ********')
    # print(zeta)
    fopt,dampopt = SSI.Cluster_Resuls(fn,zeta)
    # phi_shape = SSI.shape_modeshape(phi,fn)
    return fn,zeta ,phi,fopt,dampopt

