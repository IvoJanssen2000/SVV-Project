import math as m

"""
This program is written to convert from (z,y) coordinates on the airfoil to s-
coordinates with respect to either cell 1 or 2. The cell is chosen automatically
dependent on the coordinates inputted.
For cell 1 s = 0 is located at (z,y) = (r,0) and for cell 2 (z,y) = (0,0). In
both situations s is defined clockwise positive.

It is assumed that the airfoil is built up out of a semicircle and 2 lines, with
1 vertical spar on the intersection between semicircle and line
"""
def isWithinBounds(z,y):
    #Check if the given dimensions are within the airfoil, it is not checked
    #if the coordinates actually lie on the skin, if this is not the case
    #a wrong result will be outputted
    if((z<=r and z >= - (w-r)) and (y>=-r and y<=r)):
        return True
    else:
        return False

def calculateSCoordinate(z,y):
    if(isWithinBounds):
        #The point is in the semicircle
        if(z>0):
            theta = m.atan(y/z)
            if(y>0):
                s = r*theta
            else:
                print(s1tot)
                s = s1tot + r*theta
        #The point is in the linesection
        else:
            if(y>0):
                print("test")
                s = ymax + m.sqrt(z**2+(ymax-y)**2)
            else:
                print(-ymax-y)
                s = s2tot - h/2 - m.sqrt(z**2+(-ymax-y)**2)

        return s
    else:
        return -1


h = 0.248                    #[m]    Height of the aileron
w = 0.515                    #[m]    Cord length of the aileron
r = h/2                      #[m]    Radius of the semicircle part of the aileron
ymax = r                     #[m]    Max height of the aileron
s1tot = m.pi*r+h             #[m]    Perimeter of cell 1
slin = m.sqrt((w-r)**2+r**2) #[m]    Length of the linear part of the airfoil
s2tot = h + 2*slin           #[m]    Perimeter of cell 2
