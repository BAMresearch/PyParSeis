# -*- coding: utf-8 -*-
#!/usr/bin/python

"""
PyParseis.py

Nonlinear least squares fitting of parallel seismic first arrival data
Based on Niederleithinger 2010 (PhD thesis)

Niederleithinger BAM 8.2 2019
"""

from math import *
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from parallelseismic import *
        

# Fitfunktion
def fun(par, z, fat):
        res = PSfe(R,D,z,par[0],par[1],par[2],par[3]) - fat
        return res


# read first arrival data from ASCII text file
def load_data(filename):
        with open(filename,'r') as file:
                 dat = np.loadtxt(file)
        return dat
    
    
    
# main script
# load data
filename ='PyParseisModel-R0.3-D0.5-cp4200-cs1500-L11-off0-inc0-noise0.dat' 
#filename ='PyParseisModel-R0.3-D0.5-cp4000-cs1600-L10-off0-inc0-noise0.15.dat' 

data = load_data(filename)
z = data[:,0]
fat = data[:,1]/1000 # as fat is given in ms in file

#set parameters
R = 0.5 # pile radius in m
D = 0.3 # distance pile edge - sensor in m

# inversion parameters 
par = np.ones((5)) # starting values for inversion
par[0] = 4300.0 # concrete velocity in m/s
par[1] = 1600.0 # soil velocity  in m/s
par[2] = 12.0 # pile length in m
par[3] = 0.0001 # time offset in s
par[4] = 0.0 # inclination in deg

lower = np.ones((5)) # lower bounds for inversion (use -np.inf to disable)
lower[0] = 3000.0 # concrete velocity in m/s
lower[1] = 1000.0 # soil velocity  in m/s
lower[2] = 5     # pile length in m
lower[3] = -.001 # time offset in s
lower[4] = 0.0 # inclination in deg

upper = np.ones((5)) # upper bounds for inversion (use np.inf to disable)
upper[0] = 5000.0 # concrete velocity in m/s
upper[1] = 3000.0 # soil velocity  in m/s
upper[2] = 12.0 # pile length in m
upper[3] = .001 # time offset in s
upper[4] = 0.001 # inclination in deg


# evaluate start model
tstart = PSfe(R,D,z,par[0],par[1],par[2],par[3])

# fit 
fitresult = least_squares(fun, par, bounds=(lower,upper), verbose=1, args=(z, fat))
#evaluate fit model
tfit = PSfe(R,D,z,fitresult.x[0],fitresult.x[1],fitresult.x[2],fitresult.x[3])


# plot
cm2inch=1/2.54
#left_margin = 2.*cm2inch   # cm
#right_margin = .2*cm2inch  # cm
#figure_width = 8.*cm2inch # cm
#figure_height = 12.*cm2inch # cm
#top_margin = 0.2*cm2inch    # cm
#bottom_margin = 2.*cm2inch # cm
#box_width = left_margin + figure_width + right_margin   # cm
#box_height = top_margin + figure_height + bottom_margin # cm
#fig = plt.figure(figsize=(box_width,box_height))
fig = plt.figure(figsize=(13*cm2inch,17*cm2inch))
#plt.Axes(fig,[left_margin,bottom_margin,left_margin+figure_width,bottom_margin+figure_height])
plt.plot(fat*1000, -z, 'ok', label='data')
plt.plot(tstart*1000, -z, '-r', label='start')
plt.plot(tfit*1000, -z, '-g', label='fit')
plt.xlabel("t in ms")
plt.ylabel("z below reference in m")
plt.legend()
xmin, xmax, ymin, ymax = plt.axis()
plt.text(xmin+(xmax-ymin)/100, ymin+(ymax-ymin)/100, 'cp=%5.0f, cs=%5.0f, L=%5.2f, off=%5.3f, inc=%5.2f' % tuple(fitresult.x))
plt.savefig('parseis-fit1.png', dpi = 300)
plt.show()
#plt.close()
