from math import *
from numpy import *
import matplotlib.pyplot as plt


# Python program to calculate the spacing of the stiffeners including its width

def diststiff():

    C_a = 0.515 # C_a Chord length aileron m    (Dornier Do 228: 0.515 m)
    h_a = 24.8  # h_a Aileron height cm         (Dornier Do 228: 24.8 cm) 
    n_st = 11   # n_st Number of stiffners      (Dornier Do 228: 11)

    # Converting all lengths to meters

    h_a = h_a/100 # h_a Aileron height m 
    r = h_a/2

    # Calculation of the perimeter of the aileron

    Circ = pi*h_a/2                      # Calculation of the circumference of the semi-circle

    lh_a = C_a-h_a/2                     # Calculation of the chord length of the aileron

    ld_a = (lh_a**2 + (h_a/2)**2)**(1/2) # Calculation of the diagonal length of the aileron

    Peri = Circ+2*ld_a                   # Circumference + 2x diagonal lengths

    # Distance between each stiffener 

    d_st = Peri/n_st

    # Determining Y and Z positions of each stiffener

    Zlst = zeros((11))    
    Ylst = zeros((11))
    Zcent = zeros((11))
    Ycent = zeros((11))

    theta = d_st/r
    Zlst[0] = r
    Ylst[0] = 0

    Zlst[1] = r*cos(theta)
    Ylst[1] = r*sin(theta)

    Zlst[10] = Zlst[1]
    Ylst[10] = -Ylst[1]

    l = d_st-(pi/2-theta)*r
    delta = arctan(r/(C_a-r))

    Zlst[2] = -l*cos(delta)
    Ylst[2] = -l*sin(delta)+r

    Zlst[9] = -l*cos(delta)
    Ylst[9] = -l*cos(delta)+r

    for i in range(3,6):
        Zlst[i] = Zlst[i-1] - d_st*cos(delta)
        Ylst[i] = Ylst[i-1] - d_st*sin(delta)

    mag5 = ((Zlst[5]+0.391)**2 + Ylst[5]**2)**(1/2)
    diffd_st = d_st-mag5
    Zlst[6] = diffd_st*cos(delta)-0.391
    Ylst[6] = -diffd_st*sin(delta)


    for j in range(7,10):
        Zlst[j] = Zlst[j-1] + d_st*cos(delta)
        Ylst[j] = (Ylst[j-1] - d_st*sin(delta))

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

##    for i in range(len(Zlst)):
##        print("(",Zlst[i],",", Ylst[i],")", "\t i =", i)
    

    ##zline = [-r,0.391]
    ##yline = [0,0]
    ##plt.plot(-Zlst, -Ylst, 'o')
    ##plt.plot(zline, yline)
    ##plt.show()
        
    return Zlst, Ylst, Zcent, Ycent

def Arc(r, theta):
	""" Function to compute the y and z coordinates of the LE arc 
	Input Arguments:
		r = radius of the arc (float) [m]
		theta = range of theta values (numpy array) [rad]
	Output:
		y, z coordinates
	"""
	y = r*sin(theta);
	z = r*cos(theta);
	return z, y;

def Triangel(height, length, z):
	""" Function to compute the y and z coordinates of the trangel shape 
	Input Arguments:
		h = Height of the triangle (float) [m]
		length = Length of the airfoil - radius of arc (float) [m]
	Output:
		y, z coordinates
	"""
	dydz = -height/length;
	y = height + dydz*z;
	y1 = -height - dydz*z;
	return y, y1;

Zlst,Ylst,Zcent,Ycent = diststiff()

N = 10
h = 0.248
t = 1.1e-3
length = 0.515
theta = linspace(0, pi/2, N)
z_cord = linspace(0, length - h/2, N)
z, y = Arc(h/2, theta)
y_, y_1 = Triangel(h/2, length - h/2, z_cord)

zline = [-0.124,0.391]
yline = [0,0]
zspar = [0,0]
yspar = [0.124,-0.124]

plt.plot(-z, y, color='blue')
plt.plot(-z, -y, color='blue')
plt.plot(z_cord, y_, color='blue')
plt.plot(z_cord, y_1, color='blue')
plt.plot(-Zlst, -Ylst, 'o', color='red')
plt.plot(zline, yline, color='black')
plt.plot(zspar,yspar, color='black', linestyle='dashed')
plt.plot(-Zcent, -Ycent,"o" ,color="purple")
plt.grid(True)
plt.show()
