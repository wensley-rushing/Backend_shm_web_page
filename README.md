# shm_integration
Repository containing app integration, to automatize data driven condition assessment 

# Tutorials 

In the following links, you can find explanations and tutorials: [Google Drive Folder](https://drive.google.com/drive/folders/1GrHLS41ZESwUBPPv-Dl3kRQTudeqGdhg?usp=sharing)


## hourly_analysis.py

This code executes the tasks that need to be done each hour. The tasks are:

- **db_con.connection_and_data_retrieving:**
  Returns the engine to connect to the Postgres DB and a DataFrame object of the given dates Q1-Q2.
  
- **event_detect.get_events:**
  Returns a tuple where the first position contains an array of the picks found and the second position contains an array of the times when each event took place in the given time Q1-Q2.
  
- **modal_anl.model_analysis:**
  Returns the results of the operations made between each event detected in the `get_events` function.
  
- **time_hst.Time_history:**
  Returns the results of the operations made in text files.
  
- **rainflow.store_results_in_db:**
  Executes the rainflow operations on each file created by the `Time_history` function and saves the results in the DB.

## daily_analysis.py

This code executes the tasks that need to be done at the end of each day (Australian time). The tasks are:

- **db_conn.create_connection:**
  Returns the engine to connect to the Postgres DB.
  
- **model_upt_dly.get_model_analysis_data:**
  Returns the results of the operations made in that function and creates an OpenSees object.
  
- **model_upt_dly.create_pickle_visualization:**
  Creates the pickle and text files with the visualization given the OpenSees object.
  
        
    - crack_detec.montecarlo_simulation:
        
    - risk_mgmt.calculate_risk:
