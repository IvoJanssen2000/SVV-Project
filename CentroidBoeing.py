import numpy as np
import math as m

zcentroidlst = [ 0.1025    ,  0.06557951, -0.01831074, -0.10634515, -0.19437956,
       -0.28241397, -0.37044838, -0.45848279, -0.45848279, -0.37044838,
       -0.28241397, -0.19437956, -0.10634515, -0.01831074,  0.06557951]
ycentroidlst = [ 0.        ,  0.07877549,  0.09876497,  0.08080771,  0.06285044,
        0.04489317,  0.0269359 ,  0.00897863, -0.00897863, -0.0269359 ,
       -0.04489317, -0.06285044, -0.08080771, -0.09876497, -0.07877549]



######################## Part I - parameters as in assignment #######################################
aircraft = "B737" # Write either A320, F100, CRJ700 or Do228 (bear in mind capitals); this is used for aerodynamic loading
nst = 15   
Ca = 0.605               # m
la = 2.661               # m
x1 = 0.172               # m
x2 = 1.211               # m
x3 = 2.591               # m
xa = 0.35                # m
ha = 0.205               # m
tsp = 2.8/1000           # m
tst = 1.2/1000           # m
tsk = 1.1/1000           # m
hst = 16./1000           # m
wst = 19./1000           # m
d1 = 1.154/100           # m
d3 = 1.840/100           # m
P = 97.4*1000            # N

def areas():        #### gives an array with 15 positions. Stores  Pos. 0-10 area of stringer, pos. 11 area spar, pos.12 area semi circle, pos.13+14 area skin
    a = []
    for i in range(nst):
        a += [(wst+hst)*tst]
    ahc = 0.5 * m.pi * ha * tsk  # Area of half circle
    ask = m.sqrt((ha / 2) ** 2 + (Ca - ha / 2) ** 2) * tsk  # Area of one plane skin
    asp = ha * tsp  # Area of spar
    return (np.append(a, [asp, ahc, ask, ask]))

def z_centroids():      #### gives an array with 15 positions. it stores the z-coordinates of the centroid of each component. Pos. 0-10 stringer, pos. 11 spar, pos.12 semi circle, pos.13+14 skin
    z = zcentroidlst
    z_bar_hc = ha / np.pi
    z_bar_sk = -((Ca - ha / 2) / 2)
    z_bar_sp = 0
    return np.append(z,[z_bar_sp, z_bar_hc, z_bar_sk, z_bar_sk])

def y_centroids():      #### gives an array with 15 positions. it stores the y-coordinates of the centroid of each component. Pos. 0-10 stringer, pos. 11 spar, pos.12 semi circle, pos.13+14 skin
    return [0, 0, 0.25*ha, -0.25*ha]

areas = areas()
z_centroids = z_centroids()

def get_total_z_centroid():     #### gives z-coordinate of centroid of whole strucutre
    return sum(z_centroids*areas) / sum(areas)

zcent = get_total_z_centroid()
print(zcent)

#def get_total_y_centroid():     #### gives y-coordinate of centroid of whole strucutre
#    #return sum(y_centroids()*areas()) / sum(areas())
#    return 0
#
##offset_centroid_crosspoint_stringer = (hst/2)*hst*tst/areas()[0]
#
#z_centroids = z_centroids()
#y_centroids = y_centroids()
#
#print(z_centroids)
#
##print(z_centroids())
#
##print(get_total_z_centroid())
#
##print(diststiff()[2])
##print(diststiff()[3])
#
#print(get_total_z_centroid())
#print(get_total_y_centroid())