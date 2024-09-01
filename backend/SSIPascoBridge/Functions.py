#-------------------------------- Functions ----------------------------------#

#-- Data Analysis --#
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from scipy import signal

#-- Opensees --#

import numpy as np
import opensees.openseespy as ops
import warnings
import openseespy.postprocessing.ops_vis as opsv

import psycopg2
import sys

import pandas as pd

warnings.filterwarnings("ignore")



#-------------------------------- 1. Retrieve Info data base -----------------#


def retrieve_Data(Q1,Q2):
        try:
            cnx = mysql.connector.connect(user='python_data_source',password='minckadata2020',host='128.199.253.233',
                                        database='shm_data') 
        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print("Something is wrong with your user name or password")
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print("Database does not exist")
            else:
                print(err)
        
        
        cursor = cnx.cursor()
        

        
        query = f"SELECT * FROM `100hz_python_SIM` WHERE timestamp BETWEEN '{Q1}' and '{Q2}'"
        
        if cursor.execute(query): print("Data Saved. \n")
        
        
        data = cursor.fetchall()
        Acc = np.array([])
        A1,A2,A3,A4,A8,A5,A7 = [],[],[],[],[],[],[]
        
        for row in data:
            print("Row #: "+str(row[0]) + " AI1: "+str(row[1])+"\n")
            A1.append(row[0])
            A2.append(row[1])
            A3.append(row[2])
            A4.append(row[3])
            A5.append(row[4])
            A7.append(row[6])
            A8.append(row[7])
            
        
            
        
            
        
        
        cnx.close()  
        Acc = np.column_stack((A1,A2,A3,A4,A5,A7,A8))
        Nc = Acc.shape[1]
        for i in range(Nc):
            Acc[:,i] = signal.detrend(Acc[:,i])
        return Nc,Acc
    
    


#-------------------------------- 2. year month Date time increase -----------#

def time_Dt2_year(year,month,day,hour,minute,seconds,dt_H,dt_min):
    
    # Just use January the first, 2000
    if int(hour)>23:
        d1 = datetime(int(year), int(month), int(day)+1, int(hour)-23, int(minute), int(seconds))
    else:
        d1 = datetime(int(year), int(month), int(day), int(hour), int(minute), int(seconds))

    d2 = d1 + timedelta(hours=int(dt_H), minutes=int( dt_min))
    
    d2= d2.strftime("%Y-%m-%d %H:%M:%S")
    d1= d1.strftime("%Y-%m-%d %H:%M:%S")
    return d2,d1




#-------------------------------- 3 Date time increase ---------------------#
def Date_time_dt(current_time,duration):


    
    parsed_time = datetime.strptime(current_time, "%I:%M:%S")
    parsed_duration = datetime.strptime(duration, "%I:%M:%S")
    
    then = parsed_time + timedelta(hours=parsed_duration.hour,
                                   minutes=parsed_duration.minute,
                                   seconds=parsed_duration.second)
    
    result = then.strftime("%I:%M:%S")
    return result

#-------------------------------- 4 plot time -------------------------------#

def plt_acc_psd(Acc,fs,d2):
    Nc = Acc.shape[1]
    for i in range(Nc):
        Acc[:,i] = signal.detrend(Acc[:,i])
    
    Nfft = len(Acc[:,1]/2)+1
    plt.figure()
    plt.plot(Acc)
    plt.title(d2)
    plt.xlabel('Data Points')
    plt.ylabel('Acceleration Records')
    
    # plt.figure()
    # for i in range(Nc):
    #     plt.psd(Acc[:,i],Nfft,fs)
    
    plt.figure()
    for i in range(Nc):
        f,pxx = signal.welch(Acc[:,i],fs,nfft = Nfft,nperseg=1024)
        plt.xlim([1,int(fs/2)-1])
        plt.plot(f, 10*np.log(pxx))
        plt.title(d2)
        plt.xlabel('Data Pints')
        plt.ylabel('Power Spectral Desnsity')



#-------------------------------- 5  Save image ------------------------------#

def saveimage(fn):
    fn = fn.replace(':','_')
    root = 'Results/'
    ext = '.png'
    name = root+fn+ext
    plt.title(name)
    plt.savefig(name)
    
#-------------------------------- 6  Opensees --------------------------------#
    
    
def Op_Pasco_Bridge(A, Iy, Iz, J, E):
    ops.wipe()
    # Bridge Geometry
    L, H, z, z2 = 1000,200,175,175*2
    Nodes,conect,Supp = MG.Model_Geometry(L,H,z,z2)
    warnings.filterwarnings("ignore")
    ops.wipe()
    gama =10.5*10**-6
    g =  9800
    G =852
    MC = 0
   
       
    ops.model('basic','-ndm',3,'-ndf',6)
    z = [0,0,1]
       
    masa = np.zeros((len(Nodes),2))
        # Nodes
    for i in range(len(Nodes)):

       
        ops.node(Nodes[i,0],Nodes[i,1],Nodes[i,2],Nodes[i,3])
        masa[i, 0] = i+1
       
       
       
        ## Elementos tipo columnas
    for i in range(len(conect[:,1])):
        XYZI = ops.nodeCoord(int(conect[i,0]))
        XYZJ = ops.nodeCoord(int(conect[i,1]))
        xaxis = np.subtract(XYZJ,XYZI)  
        vecxz = np.cross(xaxis,z)
        L =np.linalg.norm(xaxis)
            #Mass
        m =A*L*gama/2/g
        MC = MC+m*2
           
        if np.linalg.norm(vecxz) == 0:
            ops.geomTransf('Linear', i,0,-1,0)
            ops.element('elasticBeamColumn', i,int(conect[i,0]),int(conect[i,1]), A, E, G, J, Iy, Iz, i)
            if  ops.nodeCoord(int(conect[i,0]))[2]>ops.nodeCoord(int(conect[i,1]))[2]:
              index = np.where(masa[:,0] == conect[i,0])[0][0]
              masa[index,1] += 2*m
            else:
              index = np.where(masa[:,0] == conect[i,1])[0][0]
              masa[index,1] += 2*m
        else:
            ops.geomTransf('Linear', i,*vecxz)
            ops.element('elasticBeamColumn', i,int(conect[i,0]),int(conect[i,1]), A, E, G, J, Iy, Iz, i)
            index = np.where(masa[:,0] == conect[i,0])[0][0]
            masa[index,1] += m
            index = np.where(masa[:,0] == conect[i,1])[0][0]
            masa[index,1] += m
   
    for i in range(len(masa)):
        ops.mass(masa[i, 0], masa[i, 1], masa[i, 1], masa[i, 1], 0, 0, 0)
    for i in range(len(Supp)):
         ops.fix(int(Supp[i]),1,1,1,1,1,1)
     
   
   
    warnings.filterwarnings("ignore")
    
    ops.timeSeries('Linear', 1)
    ops.system('BandSPD')
    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1)
    ops.analysis('Static')
    ops.analyze(1)    
   
    freq = ops.eigen('-genBandArpack', 10)
   
   
    for i in range(len(freq)):
        freq[i]  = freq[i]**0.5/(2*np.pi)
    return freq[0:2]


#-------------------------------- 7 Date time increase ---------------------#


def retrieve_Data2(Q1,Q2):
        try:
            con = psycopg2.connect(
            host="db-shm-mincka.ckensqtixcpt.ap-southeast-2.rds.amazonaws.com",
            database='shm_hailcreek',
            user='db_admin',
            password='minckA.2020')

        except psycopg2.DatabaseError as e:

            print(f'Error {e}')
            sys.exit(1)
        
        cur = con.cursor()
           
        query = f'''SELECT * FROM "100hz_python_SIM" sim WHERE sim.sampledatetime >='{Q1}' AND sim.sampledatetime < '{Q2}' ;'''
        Datao = pd.read_sql(query,con).to_numpy()
        Data =  np.delete(Datao ,[0,6,7,8,9,10,11,12,13,14,15,16],1)
        Nc= Data.shape[1]
        for i in range(Nc):
            Data[:,i] = signal.detrend(Data[:,i])
        print(Data.shape)                     
        return Nc,Data
        
        
        
        
        # query = f"SELECT *FROM `100hz_python_SIM` WHERE sim.sampledatetime >= \'2020-12-24 10:00:00\'\n AND sim.sampledatetime < \'2020-12-24 10:10:00\'\n"
        # get_ipython().run_cell_magic('time', '', 'df = pd.read_sql(\n    \'\'\'\n    SELECT *\n    FROM "100hz_python_SIM" sim\n    WHERE sim.sampledatetime >= \'2020-12-24 10:00:00\'\n        AND sim.sampledatetime < \'2020-12-24 10:10:00\'\n    ;\n    \'\'\',\n    con\n)')


        
        # query = f"SELECT * FROM `100hz_python_SIM` WHERE timestamp BETWEEN '{Q1}' and '{Q2}'"
        
        # if cur.execute(query): print("Data Saved. \n")
        
        
        # data = cur.fetchall()
        # N = len(data)
        # np.zeros((7,N))
        # Acc = np.array([])
        # A1,A2,A3,A4,A8,A5,A7 = [],[],[],[],[],[],[]
        
        # breakpoint()
        # for row in data:
        #     # print("Row #: "+str(row[0]) + " AI1: "+str(row[1])+"\n")
        #     A1.append(row[1])
        #     A2.append(row[2])
        #     A3.append(row[3])
        #     A4.append(row[4])
        #     A5.append(row[6])
        #     A7.append(row[7])
        #     A8.append(row[8])
        return Data
  
        



























    
    
