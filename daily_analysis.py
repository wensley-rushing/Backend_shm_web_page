'''
This code executes the codes that need to be done at the end of each day (Australian time). 
The codes are:
    - db_conn.create_connection:
        Returns de engine to connect to the Postgres DB.
    - model_upt_dly.get_model_analysis_data:
        Returns the results of the operations made in that function and creates and 
        opensees object.
    - model_upt_dly.create_pickle_visualization:
        Create the pkls and txt files with the visualization given the opensees object.
    - design_ratio.create_files:
        
    - crack_detec.montecarlo_simulation:
        
    - risk_mgmt.calculate_risk:
    
'''
from datetime import datetime, timedelta

import backend.design_ratio as design_ratio
import backend.crack_detection as crack_detec
import backend.database_connection as db_conn
import backend.model_updating_daily as model_upt_dly
import backend.risk_management as risk_mgmt

startTime = datetime.now()
print(startTime)
engine = db_conn.create_connection()


##Model updating
print('Doing Model Update')
fre,dampopt,phi = model_upt_dly.get_model_analysis_data (engine)
opt = model_upt_dly.opti(fre,phi)

model_upt_dly.create_pickle_visualization()
print('Done Model Update')

## Design Ratio
design_ratio.create_files()

## Crack Detection
crack_detec.montecarlo_simulation(engine)

risk_mgmt.calculate_risk()

print(datetime.now() - startTime)