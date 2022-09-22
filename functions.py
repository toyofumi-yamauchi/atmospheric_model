def load_NRLMSIS_data(filename):
    '''
    Converting the .txt file into the array
    .txt file is downloaded from https://kauai.ccmc.gsfc.nasa.gov/instantrun/msis

    '''
    # Import libraries
    import numpy as np
    import re

    # Importing the .txt file
    with open(filename,'r') as f:
        lines = f.readlines()
    print(lines[0])
    #print(lines[1])

    # Converting the lines into rows
    num_rows = np.shape(lines)[0]-1
    num_cols = np.shape(re.split(' +',lines[1]))[0]
    rows = np.zeros((num_rows,num_cols))
    for i in range(0,np.shape(lines)[0]-1):
        rows[i,:] = re.split(' +',lines[i+1])
        #print(rows[i,:])

    data = rows

    
    return data

