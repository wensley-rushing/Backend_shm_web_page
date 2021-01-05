#                               SSI COV 
#                              front end

import SSI_Backend as SSI 
import matplotlib.pyplot as plt 
import numpy as np
fn = 'INPUTS/test_1Mass.txt'
delim = ','
fs = 1/0.0667

breakpoint()
# --------------------------- 1. Load Data -----------------------------------#
Acc,Nc,N = SSI.load_data(fn,delim)
# --------------------------- 2. pot_Acc -------------------------------------#
# SSI.plot_ACC(Acc,0,Nc)
# --------------------------- 3. NexT ----------------------------------------#
IRF = SSI.NexT(Acc,fs,20,Nc)
# --------------------------- 4. blockToeplitz -------------------------------#
[U,S,V,T] = SSI.blockToeplitz(IRF)
# # --------------------------- 5. Modal id ------------------------------------#

plt.figure()
pxx,freqss = plt.psd(Acc[:,1],int(len(Acc[:,1])/2)+1,fs,color = 'blue')

fig, ax1 = plt.subplots()
color = 'tab:blue'

ax1.set_xlabel('frequency [Hz]')
ax1.set_ylabel('Power Spectral Density', color=color)

ax1.plot(freqss, 10*np.log(pxx), color=color)
plt.xlim([0,1.5])
ax1.tick_params(axis='y', labelcolor=color)
ax2 = ax1.twinx() 
color = 'tab:red'
ax2.set_ylabel('order', color=color)  # we already handled the x-label with ax1

Nmax = 22
Nmin = 4
kk=0

print('******* Identified poles **********')
fn2,zeta2,phi2,MAC,stablity_status = {},{},{},{},{}
for i in range(Nmax,Nmin,-1):
    if kk == 0:
        fn0,zeta0,phi0 = SSI.modalID(U,S,i,Nc,fs)

    else:
        fn1,zeta1,phi1 = SSI.modalID(U,S,i,Nc,fs)
        
        print(fn1)
        ax2.plot(fn1,i*np.ones(len(fn1)),'x',color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        
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
stablity_status = SSI.flip_dic(stablity_status)


fnS,zetaS,phiS,MACS = SSI.getStablePoles(fn2,zeta2,phi2,MAC,stablity_status)
fn,zeta,phi = SSI.ClusterFun(fnS,zetaS,phiS)
print('********************')
print('******** Stable Frequencies *********')
print(fnS)
print('********Stable Dampings ********')
print(zetaS)


print('********************')
print('******** cluster :D *********')
print(fn)
print('********  cluster :D  ********')
print(zetaS)




# SSI.plot_stabilization_diagram_2(fn0)
    

