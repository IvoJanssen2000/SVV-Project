import math as m

"""
This program is written to convert from (z,y) coordinates on the airfoil to s-
coordinates with respect to either cell 1 or 2. The cell is chosen automatically
dependent on the coordinates inputted.
For cell 1 s = 0 is located at (z,y) = (r,0) and for cell 2 (z,y) = (0,0). In
both situations s is defined clockwise positive.

It is assumed that the airfoil is built up out of a semicircle and 2 lines, with
1 vertical spar on the intersection between semicircle and line. Constants are
defined on the bottem of the document.
"""


def isWithinBounds(z,y):
    #Check if the given dimensions are within the airfoil, it is not checked
    #if the coordinates actually lie on the skin, if this is not the case
    #a wrong result will be outputted
    if((z<=r and z >= - (w-r)) and (y>=-r and y<=r)):
        return True
    else:
        return False


#Convert the (z,y) coordinate of a stringer to the s coordinate of either cell 1
#or cell 2.
#INPUTS:
#z [m]: the z-coordinate of the stringer
#y [m]: the y-coordinate of the stringer
#OUTPUTS:
#s [m]: the s-coordinate of the stringer, automatically taken with respect
#to the correct cell
#n [-]: The cell the stringer is in, either 1 (semicircle) or 2 (linear section)
#NOTES:
#The case z = 0 gives an unambigious result for which cell should be chosen and
#is therefore considered an invalid input. Using z = 0 can given unintended
#result. Since no stringer is placed at the spar (z=0) this should give no
#limitations to the program
def calculateSCoordinate(z,y):
    #Check if the stringer is within the airfoil, it is NOT checked if it is
    #placed on the skin, input coordinates not on the skin are likely to produce
    #wrong results.
    if(isWithinBounds(z,y)):
        #The point is in the semicircle, s is defined with respect to cell 1
        if(z>0):
            n = 1
            theta = m.atan2(y,z)
            if(y>=0):
                s = r*theta
            else:
                s = smax + r*theta
        #The point is in the linesection, s is defined with respect to cell 2
        else:
            n = 2
            if(y>=0):
                s = pi*r/2 + m.sqrt(z**2+(ymax-y)**2)
            else:
                s = smax - pi*r/2 - m.sqrt(z**2+(-ymax-y)**2)
    
        return s,n
    #If the input coordinates are not on the airfoil the value of -1,-1 will be
    #returned
    else:
        return -1,-1
#VERIFICATION RESULTS:
#Verification is done by testing one point within every different quadrant,
#as well as testing the boundary cases for
#y = 0 and testcases to test for invalid
#results, as described before it is assumed
#that z = 0 is not inputted and therefore this is not tested either.
#(z,y)=(2,1)                    returns -1                  ,-1 CORRECT
#(z,y)=(2,0)                    returns -1                  ,-1 CORRECT
#(z,y)=(-2,0)                   returns -1                  ,-1 CORRECT
#(z,y)=(0.1,2)                  returns -1                  ,-1 CORRECT
#(z,y)=(0.1,-2)                 returns -1                  ,-1 CORRECT
#(z,y)=(0.124,0)                returns 0                   ,1  CORRECT
#(z,y)=(0.062,0.10738715)       returns 0.12985249631374293 ,1  CORRECT
#(z,y)=(0.08768124,-0.08768124) returns 1.1125509542445087  ,1  CORRECT
#(z,y)=(-0.391,0)               returns 0.6049701632528961  ,2  CORRECT
#(z,y)=(-0.1955,0.062)          returns 0.3998744538877316  ,2  CORRECT
#(z,y)=(-0.1955,-0.062)         returns 0.8100658726180607  ,2  CORRECT



#CONSTANTS
pi = m.pi
h = 0.248                    #[m]    Height of the aileron
w = 0.515                    #[m]    Cord length of the aileron
r = h/2                      #[m]    Radius of the semicircle part of the aileron
ymax = r                     #[m]    Max height of the aileron
s1tot = m.pi*r+h             #[m]    Perimeter of cell 1
slin = m.sqrt((w-r)**2+r**2) #[m]    Length of the linear part of the airfoil
s2tot = h + 2*slin           #[m]    Perimeter of cell 2
smax = pi*r + 2 *slin
