import numpy as np
"""
This program is intented to calculate the torsional stiffness of an aileron.
It is assumed that the aileron geometry is symmetric about the horizontal
and made up out of a semicircle and 2 linear sections, with a spar placed
at the intersection of the semicircle and linear sections.
The calculation is made by means of compatibility and shear flow distributions.
All the constants are defined at the bottem of this program
"""


def getTorsionalStiffness():
    C1 = ((pi*r)/tSkin+h/tSpar)/(2*A1)+h/(2*A2*tSpar)
    C2 = -1*((1/(2*A2))*(h/tSpar+(2*slin)/tSkin)+h/(2*A1*tSpar))
    C3 = -1*((h/tSpar+(pi*r)/tSkin)/(2*A1))
    C4 = h/(2*A1*tSpar)

    matA = np.matrix([[2*A1,2*A2,0],
                      [C1,C2,0],
                      [C3,C4,1]])
    vec = np.matrix([[T],
                     [0],
                     [0]])

    result = np.linalg.inv(matA)*vec
    #print(result)
    J = T/(result[2])
    print(result[2][0])
    print(J)
    print(2*A1*result[0]+2*A2*result[1])
    return J




#CONSTANTS#
pi = np.pi                  #[-]
r = 0.124                   #[m]
h = 2*r                     #[m]
c = 0.515                   #[m]
w = c-r                     #[m]
slin = np.sqrt(w*w+r*r)     #[m]

A1 = (pi*r*r)/2             #[m^2]
A2 = (w*h)/2                #[m^2]

tSkin = 1.1e-3              #[m]
tSpar = 2.2e-3              #[m]

G = 28e9                    #[Pa]
T = 1                       #[Nm] (unit torque since J is a crossectional
                            #property and thuss independent of T

getTorsionalStiffness()
