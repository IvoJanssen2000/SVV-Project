import numpy as np
import math as m
from stiffener_spacing import diststiff


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


z = diststiff()[2]
y = diststiff()[3]



def areas():        #### gives an array with 15 positions. Stores  Pos. 0-10 area of stringer, pos. 11 area spar, pos.12 area semi circle, pos.13+14 area skin
    ahc = 0.5 * m.pi * ha * tsk  # Area of half circle
    ask = m.sqrt((ha / 2) ** 2 + (Ca - ha / 2) ** 2) * tsk  # Area of one plane skin
    asp = ha * tsp  # Area of spar
    ast = hst * tst + wst * tst  # Area of one stringer
    a = np.ones(nst)*ast
    return np.append(a, [asp, ahc, ask, ask])

def z_centroids():      #### gives an array with 15 positions. it stores the z-coordinates of the centroid of each component. Pos. 0-10 stringer, pos. 11 spar, pos.12 semi circle, pos.13+14 skin
    z_bar_hc = ha / m.pi
    z_bar_sk = -((Ca - ha / 2) / 2)
    z_bar_sp = 0
    a = z
    return np.append(a, [z_bar_sp, z_bar_hc, z_bar_sk, z_bar_sk])

def y_centroids():      #### gives an array with 15 positions. it stores the y-coordinates of the centroid of each component. Pos. 0-10 stringer, pos. 11 spar, pos.12 semi circle, pos.13+14 skin
    a = y
    return np.append(a, [0, 0, 0.25*ha, -0.25*ha])

def get_total_z_centroid():     #### gives z-coordinate of centroid of whole strucutre
    return sum(z_centroids()*areas()) / sum(areas())

def get_total_y_centroid():     #### gives y-coordinate of centroid of whole strucutre
    #return sum(y_centroids()*areas()) / sum(areas())
    return 0

#offset_centroid_crosspoint_stringer = (hst/2)*hst*tst/areas()[0]

print(z_centroids())
print(y_centroids())

#print(z_centroids())

#print(get_total_z_centroid())

print(diststiff()[2])
print(diststiff()[3])

print(get_total_z_centroid())
print(get_total_y_centroid())