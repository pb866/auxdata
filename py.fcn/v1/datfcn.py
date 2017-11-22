"""
datfcn
======
version 1.0

List of function focused on data processing.
contains:
    - load
    - trans
    - ReadBetween
    - order
    - xydat
"""

import sys
#______________________________________________________________________________________________

def load(file):
    """
    Function load
    =============
        
    Purpose:
        Derive matrix with sza in 0-column and j-values in folliwing columns.
    
    Variables:
    I/O:
        file:   name of input file with TUV data
        head:   indices and strings with photoreactions
        mat:    matrix with TUV ouput, first column sza in deg.,
                following columns j values in s-1
    
    internal:
        s:      counter
    
    Dependencies:
        uses:           function ReadBetween
        called from:    photMCM (main)
    """
    
    head = ReadBetween(file,"Photolysis rate coefficients, s-1","values at z =      0.500 km"," = ")
    mat  = ReadBetween(file,"sza, deg.","------------------------------------------------------------","")
    for s in range(len(head)):
        head[s][0] = int(head[s][0])
        head[s][1] = head[s][1].replace('\n', ' ').strip()
    return head, mat
#______________________________________________________________________________________________


def trans(raw):
    """
    Function trans
    ==============
    Purpose:
        Transform matrix with TUV data from strings to floats.
        
    Variables:
    I/O:
        raw:    input matrix with strings
        data:   output matrix with floats
        
    internal:
        i, j:   indices (for matrix)
    
    Dependencies:
        called from:    photMCM (main)


    to generate list of tuples instead of matrix use:
    idx = [0,1]
    for j in range(len(raw[0])):
    data.append( (raw[idx[0]][j], raw[idx[1]][j]) )
    """

    data = [[float(raw[i][j]) for j in range(len(raw[0]))] for i in range(len(raw))]
    return data
#______________________________________________________________________________________________


def ReadBetween(file,start,end,colsep):
    """
    Function ReadBetween
    ====================
        
    Purpose:
        Read input data from file between 2 different key words.
        Data is return in a matrix. Column separators can be specified for the
        split of columns.
        
    Variables:
    I/O:
        file:   file name
        start:  key word, after which (excluding the key word) data is read in
        end:    key word, before which (excluding the key word) data is read in
        colsep: column separators, use "" to leave unspecified (use whitespaces)
        outmat: output matrix with data read in between <start> and <end>
        
    internal:
        started:    switched for data between the 2 keywords
        fin:        variable for opening input file
        line:       lines read in from input file
            
        
    Dependencies:
        called from:    datfcn.load
    """

    started = False
    outmat = []
    with open(file,'r') as fin:
        for line in fin:
            if (end in line):
                started = False
            if (started):
                if (colsep != ""):
                    outmat.append(line.split(colsep))
                else:
                    outmat.append(line.split())
            if (start in line):
                started = True
    return outmat
#______________________________________________________________________________________________


def order(val):
    """
    Function order
    ==============
        
    Purpose:
        Determine the order of magnitude of a given value (val).
        Return the order as number and a factor 10^(order).
        
    Variables:
    I/O:
        val:    input value (any number)
        ord:    order of magnitude
                (negative integers for values < 0, 0 for 0, positvie values for values > 0)
        mult:   factor 10^(order)
        
    Dependencies:
        uses:           numpy
        called from:    photMCM (main)
    """
    import numpy as np
    if (val != 0):
        ord  = np.floor(np.log10(np.abs(val)))
    else:
        ord = 0
    mult = 10**ord
    return ord, mult
#______________________________________________________________________________________________

def xydat(data,rxn,co,om,fout,gp,ffcn,scen):
    """
    Function xydat
    ==============
    
    Purpose:
        Retrieve x and y data for curve fitting (threshold can be introduced to cut off values).
        
    Variables:
    I/O:
        co:             Threshold (cut-off parameter for sza) for TUV input data in deg.
        data:           matrix with TUV data
        rxn:            header of TUV matrix (used to determine size of matrix and print progress)
        fout,gp,ffcn:   identifiers for output files
        om:             matrix with all orders of magnitudes for parameter l of fitting curve
        scen:           scenario name
        
    internal:
        x,y:            counters/indices
        xdata:          data of independent variable (sza) for curve fitting
        ydata:          data of dependent variable (j values) for curve fitting
    
    Dependencies:
        uses:           numpy, fitTUV, fitStat, pfit
        called from:    photMCM (main)
    """
    # import functions
    import numpy as np
    from fitfcn import *
    from pltfcn import *
    xdata =[]
    
    # declare x data for fit (sza)
    for x in range(len(data)):
        if (data[x][0] < np.deg2rad(co)):
            xdata = np.append(xdata,data[x][0])
# loop over photoreactions
    for y in range(len(rxn)):
        print rxn[y][1] # progress output to screen
        ydata = []
    
        for x in range(len(data)):
            if (data[x][0] < np.deg2rad(co)):
                ydata = np.append(ydata,data[x][y+1])
        if (ydata.sum() == 0.):
            print "Column contains only zeros.\nSkipping data processing."
            continue
        # least square curve fitting:
        p,l,ol,m,n,cov,infodict = fitTUV(xdata,ydata)
        # calculation of statistical data:
        rsquared,ss_tot,rmse,el,em,en = fitStat(xdata,ydata,p,cov,infodict)
        # printing output:
        pfit(fout,gp,ffcn,rxn,y,p,om,l,m,n,ol,el,em,en,rsquared,ss_tot,rmse,scen)
    return None

#______________________________________________________________________________________________