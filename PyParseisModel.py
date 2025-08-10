# -*- coding: utf-8 -*-
#!/usr/bin/python

"""
PyParseisModel.py

Calculate synthetic parallel seismic first arrival data
Based on Niederleithinger 2010 (PhD thesis)

Niederleithinger BAM 8.2 2020
"""

from math import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as pltext

from parallelseismic import *
        
# main script

#set parameters
R = 0.3 # pile radius in m
D = 0.5 # distance pile edge - sensor in m

z0 = 1.0 # shallowest sensor position in m (positive downwards)
dz = 0.5  # sensor interval in m 
z1 = 15 # deepest sensor position in m

z=np.arange(z0,z1*1.001,dz) # z1*1.001 to include z1 in array

cp = 4200.0 # concrete velocity in m/s
cs = 1500.0 # soil velocity  in m/s
L = 11.0 # pile length in m
off = 0.000 # time offset in s
inc = 0.0 # inclination in deg

# evaluate start model
tcalc = PSfe(R,D,z,cp,cs,L,off)

#add noise 
noiseamp = 0.00/1000  # noise amplitude in s
iz = z.size
noise=np.random.random_sample((iz,))
noise= (noise-0.5)*noiseamp
tcalc =tcalc+noise

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
plt.plot(tcalc*1000, -z, 'ok', label='data')
plt.xlabel("t in ms")
plt.ylabel("z below reference in m")
plt.legend()
xmin, xmax, ymin, ymax = plt.axis()
plt.text(xmin+(xmax-ymin)/100, ymin+(ymax-ymin)/100, 'cp=%5.0f, cs=%5.0f, L=%5.2f, off=%5.3f, inc=%5.2f' % (cp,cs,L,off,inc))
plt.savefig('parseis-fit.png', dpi = 300)
plt.show()
plt.close()

# write first arrival data to ASCII text file

filename='PyParseisModel-R0.3-D0.5-cp4200-cs1500-L11-off0-inc0-noise0.dat' 
f = open(filename,'w')
for i in range (0,iz):
    f.write(str(z[i])+' '+str(tcalc[i]*1000)+'\n')
f.close()

 