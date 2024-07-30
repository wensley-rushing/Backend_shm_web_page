# shm_integration
Repository containing app integration, to automatize data driven condition assessment 

# Tutorials 

In the following links, you can find explanations and tutorials: [Google Drive Folder](https://drive.google.com/drive/folders/1GrHLS41ZESwUBPPv-Dl3kRQTudeqGdhg?usp=sharing)




## hourly_analysis.py

This code executes the codes that need to be done each hour. The codes are:
    - db_con.connection_and_data_retrieving:
        Returns de engine to connect to the Postgres DB and a DataFrame object
        of the given dates Q1-Q2.
    - event_detect.get_events:
        Returns a tuple that has in the first position an array of the picks found
        and in the second an array of the times in which each event took place in 
        the given time Q1-Q2.
    - modal_anl.model_analysis:
        Returns the results of the opetarions made between each event that was
        detected in the function get_events.
    - time_hst.Time_history:
        Returns the results of the operations made in txt files 
    - rainflow.store_results_in_db:
        Executes the rainflow operations in each file that was created in the 
        Time_history function and saves the results in the DB.
        
## daily_analysi.py

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
