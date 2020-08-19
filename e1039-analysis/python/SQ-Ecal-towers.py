# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:36:44 2020

@author: sergi
"""


from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.stats import halfnorm


import uproot
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
import os,sys
import pandas as pd
from matplotlib import rcParams



def symmetricize(arr1D):
    ID = np.arange(arr1D.size)
    return arr1D[np.abs(ID - ID[:,None])]




directory = 'New Slide Files Brem'
indexer=0
for filename in os.listdir(directory):
    ifile=filename
    name=filename.split("_")
    Mass=name[3]+" GeV"
    Type=name[2]
    truth = uproot.open("New Slide Files Brem\%s"%ifile)["Truth"]
    truthdf = truth.pandas.df(flatten=False)
    
    indexer+=1
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', 3)
    print(truthdf)
    
    
    l=truthdf["hit_edep"]
    ID=truthdf["hit_detID"]
    eID=truthdf["hit_elmID"]
    E=truthdf["ge"]
    flat_edep = []
    elmID=[]
    y_pos=[]
    x_t=[]
    y_t=[]
    new=[]
    for e in range(len(ID)):
        for j in range(len(ID[e])):
            if int(ID[e][j])==100:
                flat_edep.append(l[e][j])
                elmID.append(eID[e][j])
                new.append(eID[e][j])
        if len(new) == 0:
            flat_edep = []
            elmID=[]
            x_t=[]
            y_t=[]
            new=[]
            continue
        
        new=[]
        print("Event "+str(e)+": "+str(len(elmID))+" hits")
        z = np.array(flat_edep)
        con=(0.00051099895/0.00064903)
        z=z*con
        x_t=[]
        y_t=[]
        existID=[]
        for v in elmID:
            y=v%36
            x=v//36
            x_t.append(x)
            y_t.append(y)


        plt.hist2d(x_t, y_t, bins=72, density=False)    
        plt.xlabel("x")
        plt.ylabel('y')
        cbar=plt.colorbar()
        cbar.set_label('Counts')
        plt.title("Unweighted Counts for "+Type+": "+Mass+" \n Event "+str(e))
        plt.show()
        flat_edep = []
        elmID=[]
        x_t=[]
        y_t=[]
    