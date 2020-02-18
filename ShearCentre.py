from numpy import*
from matplotlib import pyplot as plt

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


N = 10;
h = 0.248;
t = 1.1e-3;
length = 0.515;
theta = linspace(0, pi/2, N);
z_cord = linspace(0, length - h, N);
z, y = Arc(h/2, theta);
y_, y_1 = Triangel(h/2, length - h, z_cord);

plt.plot(-z, y, 'x');
plt.plot(-z, -y, 'x');
plt.plot(z_cord, y_, 'x');
plt.plot(z_cord, y_1, 'x');
plt.grid(True);
plt.show();

