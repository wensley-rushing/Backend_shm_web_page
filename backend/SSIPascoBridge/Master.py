
import Functions as Fn
import SSI_COV_AD as SSICOV
import numpy as np
import matplotlib.pyplot as plt
F = {}
z = {}
fs =100
#-------------------------------- 2. year month Date time increase ------------#
for j in range(3,4):
    for k in range(5,40,5):
        year,month,day,hour,minute,seconds = 2020.,12.,29.,int(j),int(k),0.
        dt_H,dt_min = 0, 2.
        
        d2,d1 = Fn.time_Dt2_year(year,month,day,hour,minute,seconds,dt_H,dt_min)
        
        
        # #-------------------------------- 1. Retrieve Info data base -------#
        Nc,Acc = Fn.retrieve_Data2(d1,d2)

        # #-------------------------------- SSI COV --------------------------#


        fn,zeta ,phi,fopt,dampopt = SSICOV.SSI_COV_AD(Acc,fs,3,Nc,30,5)
        print(fopt)



                           
