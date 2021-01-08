from datetime import datetime, timedelta

import backend.design_ratio as design_ratio
import backend.crack_detection as crack_detec
import backend.database_connection as db_conn
import backend.model_updating_daily as model_upt_dly
import backend.risk_management as risk_mgmt

startTime = datetime.now()

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