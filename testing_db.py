from datetime import datetime, timedelta

import backend.database_connection as db_con


date_now = datetime.now() + timedelta(hours = 10)
year = date_now.year
month = date_now.month
day = date_now.day
hour = date_now.hour
Q2 = datetime(year, month, day, hour, 0, 0)
Q1 = Q2 - timedelta(hours = 1)


print(datetime.now())
print('Loading data')

accelerations_df, engine = db_con.connection_and_data_retrieving(Q1,Q2)

print(len(accelerations_df))