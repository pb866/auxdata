"""
fitfcn
======
version 1.0

List of function focused on the parameterisation and curve fitting.
contains:
    - phot
    - residuals
    - fitTUV
    - fitStat
"""

import sys
#______________________________________________________________________________________________

def phot(p,x):
    """
        Function phot
        =============
        
        Purpose:
            Parameterisation for photoreactions in MCM.
            Parameterise sza dependent j values from TUV model calculations.
        
        Variables:
        I/O:   
            p(l,m,n):   array of parameter l (in s-1), m, n
            x:          sza (independent variable) in parameterisation
            jval:       j values (dependent variable in parameterisation)
        
        internal:
        
        
        Dependencies:
        uses:           numpy
        called from:    photMCM (main)
        """
    l,m,n=p
    import numpy as np
    jval = l*((np.cos(x))**m)*np.exp(-n/np.cos(x))
    return jval
#______________________________________________________________________________________________


def residuals(p,xdata,jval):
    """
        Function residuals
        ==================
        
        Purpose:
            Calculate the residuals: distance | true (TUV) values - calculated (fit) values |
        
        Variables:
        I/O:
            p(l,m,n):   fit paramters l (in s-1), m, n
            xdata:      sza (values of independent variable)
            jval:       j values (values of dependent variable)
        
        Dependencies:
            uses:           fitfcn.phot
            called from:    photMCM (main), scipy.optimize.leastsq
        """
    return jval - phot(p,xdata)
#______________________________________________________________________________________________


def fitTUV(xdata,ydata):
    """
    Function fitTUV
    ===============
        
    Purpose:
        Curve of TUV data to parameterisation for photolysis in MCM.
        j = l*(cos(x))^m*exp(-n*sec(x))
        
    Variables:
    I/O:
        xdata,ydata:    x-,y-data for curve fit (sza and j values)
        p:              optimised parameters from least square fit
        cov:            covariance matrix from least square fit
        infodict:       further statistical data from least square fit
        mesg, ier:      not used (see description of leastsq function)
        
    internal:
        Param           list with parameters for curve fitting
        p_guess:        initial parameters for curve fit
        
    Dependencies:
        uses:           residuals,numpy,scipy.optimize.leastsq, collections.namedtuple
        called from:    datfcn.xydat
    """
    # load functions
    import numpy as np
    import collections #
    from scipy.optimize import leastsq # for curve fit
    from datfcn import order
    
# curve fitting
    Param=collections.namedtuple('Param','l m n')
    p_guess=Param(l=ydata[0]*np.exp(-0.8),m=0.8,n=0.25)
    p,cov,infodict,mesg,ier = leastsq(residuals,p_guess,args=(xdata,ydata),full_output=True)
    p=Param(*p)
    l = p[0]
    ol, ml = order(l)
    l = l/ml
    m = p[1]
    n = p[2]
    return p,l,ol,m,n,cov,infodict
#______________________________________________________________________________________________


def fitStat(xdata,ydata,p,cov,infodict):
    """
    Function fitStat
    ================
        
    Purpose:
        Derive statistical data from least square fit of TUV data for MCM parameterisations.
        
    Variables:
    I/O:
        xdata,ydata:    x-,y-data for curve fit (sza and j values)
        p:              optimised parameters from least square fit
        cov:            covariance matrix from least square fit
        infodict:       further statistical data from least square fit
        rsquared:       correlation coefficient
        ss_tot:         standard deviation
        rmse:           root mean square error
        el,em,en:       confidence of parameters
        
    internal:
        ss_err:         for calculation of R2 and RMSE
        dof:            degrees of freedom for calculation of RMSE
        s_sq:           for calculation of confidences of parameters
        err:            list of all confidences
        
    Dependencies:
        uses:           residuals,datfcn.order,numpy
        called from:    datfcn.xydat
        """
    # load functions
    import numpy as np
    from datfcn import order
    
    
    # R2:
    ss_err=(infodict['fvec']**2).sum()
    ss_tot=((ydata-ydata.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    
    # RMSE:
    dof=len(xdata)-len(p)
    rmse=np.sqrt(ss_err/dof)
    
    # Confidence:
    
    if (len(ydata) > len(p)) and cov is not None:
        s_sq = (residuals(p, xdata, ydata)**2).sum()/(len(ydata)-len(p))
        cov = cov * s_sq
    else:
        cov = inf
    
    err = np.sqrt(np.diag(cov))
    el = err[0] / order(p[0])[1]
    em = err[1]
    en = err[2]

    return rsquared,ss_tot,rmse,el,em,en
#______________________________________________________________________________________________