import numpy as np
import pandas as pd
import codecs

def read_deeper(input_file:str ,*,drop:bool = True):
    d = pd.read_csv(input_file,names = ['lat','lon','depth','t'])
    if drop:
        d = d[d.lat != 0]
        d = d[d.lon != 0]
        d = d[d.t != 0]
        d.index = range(len(d))
        unique_t,unique_index = np.unique(d.t,return_index=True)
        d = d.loc[unique_index]
        d = d.sort_values(by='t')
        d.index = range(len(d))
    return d


