from numpy import *
import matplotlib.pyplot as plt


# Python program to calculate the spacing of the stiffeners including its width

def diststiff(C_a, h_a, n_st):
    """Function to compute the distance between each stiffener.
    Input arguments:
        C_a = Chord length [m]
        h_a = Aileron height [m]
        n_st = number of stiffeners [-]
    Output arguments:
        d_st = Distance between each stiffener [m]"""
    
    # Calculation of the perimeter of the aileron

    Circ = pi*h_a/2                      # Calculation of the circumference of the semi-circle
    lh_a = C_a-h_a/2                     # Calculation of the chord length of the aileron
    ld_a = (lh_a**2 + (h_a/2)**2)**(1/2) # Calculation of the diagonal length of the aileron
    Peri = Circ+2*ld_a                   # Circumference + 2x diagonal lengths

    # Distance between each stiffener 

    d_st = Peri/n_st

    return d_st

def stiffcoord(d_st, r, C_a):
    """Function to compute the Z and Y coordinates of each stiffener.
    Input arguments:
        d_st = Distance between each stiffener [m]
        C_a = Chord length [m]
        r = Radius of the semi-circle [m] 
    Output arguments:
        Zlst = Numpy array of the z-coordinates of the stiffeners [m]
        Ylst = Numpy array of the y-coordinates of the stiffeners [m]
        """
    
    # Creating 11x1 array of zeroes for the Y-axis and Z-axis for the stiffeners
    
    Zlst = zeros((11))    
    Ylst = zeros((11))
    
    # Assign the first stiffener at the leading edge on the chord line

    Zlst[0] = r 
    Ylst[0] = 0
    
    # Finding the arc angle
    
    theta = d_st/r

    # Assign the coordinates of the 2nd stiffener going clockwise from the first
    
    Zlst[1] = r*cos(theta)
    Ylst[1] = r*sin(theta)

    # Assign the coordinates of the 3rd stiffener transitioning from the semi-circle to the straight line
    
    l = d_st-(pi/2-theta)*r
    delta = arctan(r/(C_a-r))

    Zlst[2] = -l*cos(delta)
    Ylst[2] = -l*sin(delta)+r

    # Assign the coordinates of the 4-6th stiffener on the top straight line
    
    for i in range(3,6):
        Zlst[i] = Zlst[i-1] - d_st*cos(delta)
        Ylst[i] = Ylst[i-1] - d_st*sin(delta)

    
    # Assign the coordinates of the 7th stiffener transitioning from the top straight line to the bottom straight line
    
    mag_5 = ((Zlst[5]+(C_a-r))**2 + Ylst[5]**2)**(1/2)  # Calculating the hypotenuse of the top straight line between the chordline and the location of the 5th stiffener
    diffd_st = d_st-mag_5                               # Calculating the hypotenuse of the bottom straight line between the chordline and the location of the 6th stiffener
    Zlst[6] = diffd_st*cos(delta)-(C_a-r)               # Calculating Z coordinate of the 6th stiffener
    Ylst[6] = -diffd_st*sin(delta)                      # Calculating the Y coordinate of the 6th stiffener

    # Assign the coordinates of the 8-10th stiffener on the bottom straight line
    
    for j in range(7,10):
        Zlst[j] = Zlst[j-1] + d_st*cos(delta)
        Ylst[j] = (Ylst[j-1] - d_st*sin(delta))

    # Assign the coordinates of the 11th stiffener since it's symmetric along the semi-circle

    Zlst[10] = Zlst[1]
    Ylst[10] = -Ylst[1]

##    for i in range(len(Zlst)):
##        print("(",Zlst[i],",", Ylst[i],")", "\t i =", i)
##
##
##    zline = [-r,0.391]
##    yline = [0,0]
##    plt.plot(-Zlst, -Ylst, 'o')
##    plt.plot(zline, yline)
##    plt.show()
        
    return Zlst, Ylst

def stiffcentcoord(Zlst, Ylst):
    """Function to compute the Z and Y coordinates of each stiffener from their centroid lcoation.
    Input arguments:
        Zlst = Numpy array of the z-coordinates of the stiffeners [m]
        Ylst = Numpy array of the y-coordinates of the stiffeners [m]
    Output arguments:
        Zcent = Numpy array of the z-coordinates of the stiffeners from centroid location [m]
        Ycent = Numpy array of the y-coordinates of the stiffeners from centroid location [m]
        """

    # Creating 11x1 array of zeroes for the Y-axis and Z-axis for the stiffeners centroid locations
    
    Zcent = zeros((11))
    Ycent = zeros((11))

    Zcent[0] = Zlst[0]-2.5/1000
    Ycent[0] = Ylst[0]

    Zcent[1] = Zlst[1]-1.58/1000
    Ycent[1] = Ylst[1]-1.94/1000

    Zcent[10] = Zlst[10]-1.58/1000
    Ycent[10] = Ylst[10]+1.94/1000

    for i in range(2,6):
        Zcent[i] = Zlst[i]+0.76/1000
        Ycent[i] = Ylst[i]-2.38/1000

    for i in range(6,10):
        Zcent[i] = Zlst[i]+0.76/1000
        Ycent[i] = Ylst[i]+2.38/1000

    return Zcent, Ycent

def Arc(r, gamma):
	""" Function to compute the y and z coordinates of the LE arc 
	Input Arguments:
		r = radius of the arc (float) [m]
		gamma = range of gamma values (numpy array) [rad]
	Output:
		y, z coordinates
	"""
	y = r*sin(gamma)
	z = r*cos(gamma)
	
	return z, y

def Triangle(height, length, z):
	""" Function to compute the y and z coordinates of the triangle shape 
	Input Arguments:
		h = Height of the triangle (float) [m]
		length = Length of the airfoil - radius of arc (float) [m]
	Output:
		y, z coordinates
	"""
	dydz = -height/length
	y = height + dydz*z
	y1 = -height - dydz*z
	
	return y, y1

# Assigning all inputs

N = 10          # Number of points
t = 1.1e-3      # Thickness of the skin [m]       (Dornier Do 228: 1.1 mm)
C_a = 0.515     # Chord length aileron [m]        (Dornier Do 228: 0.515 m)
h_a = 24.8/100  # Aileron height [m]              (Dornier Do 228: 24.8 cm) 
n_st = 11       # Number of stiffners             (Dornier Do 228: 11)
r = h_a/2       # Radius of the semi-circle [m] 


gamma = linspace(0, pi/2, N)            # Numpy array of all gamma values used for the leading edge arc
z_coord = linspace(0, C_a - h_a/2, N)   # Numpy array of all Z coordinates for the triangle

# Assigning all outputs

z, y = Arc(h_a/2, gamma)

y_, y_1 = Triangle(h_a/2, C_a - h_a/2, z_coord)

d_st = diststiff(C_a, h_a, n_st)

Zlst, Ylst = stiffcoord(d_st, r, C_a)

Zcent, Ycent = stiffcentcoord(Zlst, Ylst)

# Creating lines for plotting the chord line & the spar

zline = [-0.124,0.391]
yline = [0,0]
zspar = [0,0]
yspar = [0.124,-0.124]

# Plotting the aileron geometry

plt.plot(-z, y, color='blue')        # Plotting the first quadrant of the semi-circle
plt.plot(-z, -y, color='blue')       # Plotting the second quadrant of the semi-circle
plt.plot(z_coord, y_, color='blue')  # Plotting the top straight line along the chord
plt.plot(z_coord, y_1, color='blue') # Plotting the bottom straight line along the chord 

# Plotting the stiffeners 

plt.plot(-Zlst, -Ylst, 'o', color='red')
plt.plot(-Zcent, -Ycent,"o" ,color="purple")

# Plotting the chord line & the spar

plt.plot(zline, yline, color='black', linestyle='dashed')
plt.plot(zspar,yspar, color='black')

plt.grid(True)
plt.show()
