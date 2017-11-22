"""
    List of function focused on plotting and writing file output.
    contains:
    - gpscat
    - pfit
"""

import sys
#______________________________________________________________________________________________
    
def gpscat(gpfile,rxn,data):
    """
        Function gnuinp
        ===============
        
        Purpose:
            Write gnuplot input file with data for scatter plots of TUV data.
        
        Variables:
        I/O:
            gpfile:     gnuplot file with data for scatter plots
            rxn:        matrix with indices and labels of available photoreactions
            data:       matrix with sza-dependent j values from TUV model
            om:         matrix with all orders of magnitudes (for sza column 0 dummy inserted)
        
        internal:
            h,x,y:      counters/indices
            o1,o2:      order of magnitude of maximum value in each j data column
        
        Dependencies:
        uses:           numpy, datfcn.order
        called from:    photMCM (main)
        """

    # initilise data
    import numpy as np
    from datfcn import order
    om    = []

    # open text file
    with open(gpfile,'w+') as gnu:
        # write header
        gnu.write("# sza / rad ")
        for h in range(len(rxn)):
            gnu.write("\t%s" % (rxn[h][1]))
        gnu.write("\n")
        # determine order of magnitude from maximum (first) value
        # for nicer output
        for y in range(len(rxn)):
            o1, o2 = order(data[0][y+1])
            om = np.append(om,o2)
        om = np.insert(om,0,1.)
            # write data matrix
        for x in range(len(data)):
            for y in range(len(rxn)+1):
                if (y == 0):
                    data[x][0] = np.deg2rad(data[x][0])
                gnu.write("\t%.3f" % (data[x][y]/om[y]))
            gnu.write("\n")
    return rxn, data, om
#______________________________________________________________________________________________

def pfit(fout,gp,ffcn,rxn,y,p,om,l,m,n,ol,el,em,en,rsquared,ss_tot,rmse,scen):
    """
    Function pfit
    =============
        
    Purpose:
        Write gnuplot input file with data for scatter plots of TUV data.
        
    Variables:
    I/O:
        fout,gp,ffcn:   identifiers for output files
        rxn:            matrix with indices and labels of available photoreactions
        p:              matrix with all optimised paramters from least square fit
        l,m,n:          optimised paramters from least square fit
        el,em,en:       confidence of parameters
        ol:             order of magnitude of parameter l
        om:             matrix with all orders of magnitudes (for sza column 0 dummy inserted)
        rsquared:       correlation coefficient
        ss_tot:         standard deviation
        rmse:           root mean square error
        y:              index (for addressing positions in rxn-matrix)
        now:            variable for present time
        scen:           scenario name
        
    Dependencies:
        uses:           datetime
        called from:    datfcn.xydat
    """
    
    import numpy as np
    from datetime import datetime as dt
    now = dt.now()
    # write data in table in file <scen>.dat
    fout.write("(%.3f%s%.3f)E%i\t %.3f%s%.3f \t %.3f%s%.3f \t  %.3e \t %.4f \t%s\n" %(l,u"\u00B1",el,ol,m,u"\u00B1",em,n,u"\u00B1",en,rmse,rsquared,rxn[y][1]))
    # write gnuplot commands into file gnu.plt
    gp.write("j%i(x) = %.3f*(cos(x))**(%.3f)*exp(-(%.3f)/cos(x))\n" % (rxn[y][0],p[0]/om[y+1],m,n))
    gp.write("set ylabel \"j(%s) / 10^{%i} s^{-1}\"\n" % (rxn[y][1], int(np.log10(om[y+1]))))
    gp.write("set y2label \"created %s.%s.%s, %s:%s\" font \",16\"\n" % (now.day,now.month,now.year,now.hour,now.minute))
    gp.write("plot \'gnu.dat\' u 1:%i ti 'TUV output' w p ls 1,\\\n     j%i(x) ti \'Least square fit\' w l ls 2\n\n" % (rxn[y][0]+1,rxn[y][0]))
    # write gnuplot functions into separate file, which is saved for later processing
    ffcn.write("j%i_%s(x) = %.3f*(cos(x))**(%.3f)*exp(-(%.3f)/cos(x)) # order of magnitude: %i\n"
             % (rxn[y][0],scen,p[0]/om[y+1],m,n,int(np.log10(om[y+1]))))

    return None

#______________________________________________________________________________________________

