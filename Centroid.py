import numpy as np
import math as m

######################## Part I - parameters as in assignment #######################################
aircraft = "Do228" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
Ca = 0.515       # m
la = 2.691       # m
x1 = 0.174       # m
x2 = 1.051       # m
x3 = 2.512       # m
xa = 0.30       # m
ha = 0.248       # m
tsk = 1.1/1000  # m
tsp = 2.2/1000  # m
tst = 1.2/1000  # m
hst = 1.5/100  # m
wst = 3.0/100  # m
nst = 11      # -
d1 = 1.034/100       # m
d3 = 2.066/100       # m
theta = m.radians(25)  # rad
P = 20.6*1000  # N


################ Calculate all areas #########


Ahc = 0.5*m.pi*ha*tsk # Area of half circle
Ask = m.sqrt((ha/2)**2 + (Ca - ha/2)**2) * tsk # Area of one plane skin
Asp = ha*tsp # Area of spar
Ast = hst*tst + wst*tst # Area of one stringer

Atot = Ahc + 2*Ask + Asp + nst*Ast

############ Calculate respective z-value of centroid of each component ############

#z_bar_hc = (4 * ha/2) / 3*m.pi  #still wrong
z_bar_sk = -((Ca-ha/2)/2)
z_bar_sp = 0
offset_centroid_crosspoint_stringer = (hst/2)*hst*tst/(Ast)


########### Calculate total z-value of centroid ############

#z_bar = (z) / Atot

print(offset_centroid_crosspoint_stringer)