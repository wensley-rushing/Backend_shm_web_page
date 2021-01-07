import pandas as pd
import numpy as np
from scipy import stats, signal, sparse

def init_limits(df):
    _dict = dict({key:None for key in df.columns[1:]})
    return _dict


def update_limits(df, dict_limits, cols = None,index = True,n_bins = 100,ci = 4):
    '''
    Compute the statistics mean, variance and standard deviation of a given dataframe columnwise.

    Parameters
    ----------
    df : Dataframe
        Input data.
        
    cols : list of str, optional
        list with the name of the columns in the dataframe for wich calculates the statistics, by default
        it takes all the columns except of the first one if index is set False.
        
    index : bool, optional
        If set True it ignores the first column when calculating the statistics, default True.
        
    n_bins : int, optional
        It defines the number of equal-width bins in the given range (100, by default).

    '''
    if not index:
        starts = 0
    else:
        starts = 1
        
    if cols is None:
        cols = df.columns[starts:].to_list()
        
    for col in cols:
        _stats =  stats.rv_histogram(np.histogram(df[col].values,bins = n_bins)).stats()
        mu,sd = _stats[0].item(),np.sqrt(_stats[1].item())
        new = mu + ci*sd
        if dict_limits[col] is None:
            dict_limits[col] = new
        else:
            old = dict_limits[col]
            dict_limits[col] = old + 0.9*(new-old)
    
    return dict_limits

def get_events(df,delta = 80,time_unit = 's'):
    '''
    Find the events of interest in the given data.
    
    Parameters
    ----------
    df : Dataframe
        Input data
    delta : int, optional
        
    time_unit: {'s'}
    
    Retunrs
    -------
    Peaks
    
    events_date
    
    '''
    nrows,ncol = df.shape
    peaks = sparse.lil_matrix((nrows,ncol-1))
    #peaks = np.zeros((nrows,ncol-1))
    isEvent = np.zeros(nrows,dtype=bool)
    
    dict_limits = init_limits(df)
    dict_limits = update_limits(df,dict_limits)
    
    for col,sensor in enumerate(df.columns[1:]):
        _max = dict_limits[sensor]
        rows,vals = signal.find_peaks(df[sensor].values,_max)
        peaks[rows,col] = vals['peak_heights']
        isEvent[rows] = True
        
    times = df.loc[isEvent,'sampledatetime'].values
    time_events = list()
    
    i = 0
    while i<times.size:
        j = 0
        while i+j+1<times.size and times[i+j+1]-times[i+j]<np.timedelta64(delta,'s'):
            j+=1
        time_events.append((times[i],times[i+j]))
        i += j+1
    return peaks.tocsc(),np.array(time_events)