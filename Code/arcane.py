#this file contains arcane functions that should not exist
#author: Nicolaas van der Voort
#July 17, 2020
#HHU, Dusseldorf


import pandas as pd

def renameDataFrameForMargarita(DataFrame):
    """takes a DataFrame and alters the column names to be Margarita readable"""
    renames = {
        'NG' : 'Green Count Rate (KHz)',
        'NR' : 'Red Count Rate (KHz)',
        'NY' : 'Yellow Count Rate (KHz)',
        'posxG' : 'peak_x_cI',
        'posyG' : 'peak_y_cI',
        'posxY' : 'peak_x_cII',
        'posyY' : 'peak_y_cII'
    }
    for oldname, newname in renames.items():
        try:
            DataFrame.rename(columns = {oldname : newname}, inplace = True)
        except:
            print('name %s not found' %oldname)
    return DataFrame