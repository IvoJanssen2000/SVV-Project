from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import inv, solve, eigvals, det, norm
import time

def Coordinate(N, l_ref, dir):
	""" Function to compute the coordinate values along the Span and Chord wise direction
	Input Arguments:
		N = Number of nodes in the given direction (int)
		l_ref = Physical Length in the given direction (float)
		dir = Coordinate direction, "x" or "z" (string)
	Output:
		pos = 1D array containing nodes in the direction (numpy array)
	"""
	if dir == "x":
		c = 1;
	elif dir == "z":
		c = -1;
	pos = zeros((1, N))[0];
	for i in range(1, N + 2):
		theta_i = (i - 1)*pi/N;
		theta_j = i*pi/N;
		if i != N + 1:
			#print(i, N);
			pos[i - 1] = c*0.5*(l_ref/2*(1 - cos(theta_i)) + l_ref/2*(1 - cos(theta_j)));
	return pos;

def Cosine_Sampler(a, b, n):
	""" Function to sample n points in [a, b] with non-uniform (Cosine) distribution
	Input Arguments:
		a = Start point (float)
		b = End point (float)
		n = Number of samples (int)
	Output:
		xi = 1D array containing n points in [a, b]
	"""
	ksi = zeros((1, n));
	for i in range(1, n + 1):
		ksi[0][i - 1] = cos((2*i - 1)*pi/(2*n));
	xi = (a + b)/2 + (b - a)/2*ksi;
	return xi[0];

## Quadrature
def GenMatrix(N, x):
	""" Function to Generate Left Matrix (A) for quadrature
	Input Arguments:
		N = Degree of Precision (int)
		x = Quadrature nodes (numpy float array)
	Output:
		A = Inverse of the left side matrix
	"""
	A = ones((N, N));
	for i in range(1, N):
		A[i, :] = x**i;
	return inv(A);

def Quadrature_weights(x):
	""" Function to compute the weights for the Quadrature rule
	Input Arguments:
		x = Quadrature nodes (numpy float array)
	Output:
		w = Weights of the Quadrature rule
	"""
	#x = Cosine_Sampler(a, b, N);
	N = len(x);
	a = x[0]; b = x[-1];
	A_i = GenMatrix(N, x);
	B = zeros((N, 1));
	for i in range(1, N + 1):
		B[i - 1] = (b**i - a**i)/i;
	w = A_i.dot(B);
	return w.T[0];

def integrate(f, wx, wy):
	""" Function to Compute 2D integral using Quadrature rule
	Input Arguments:
		f = Matrix containing function values at the nodes (numpy array)
		wx = Quadrature weights in x (numpy array)
		wy = Quadrature weights in y (numpy array)
	Output:
		I = Value of Integral
	"""
	I = f.dot(wx).dot(wy);
	return I;

## Radial Basis Function Interpolation
#def RBF(x1, x2):											# Old Function, Took too long
#	""" Function to define Radial Basis Function
#	Input Arguments:
#		x1 = node 1 (float or numpy array)
#		x2 = node 2 (float or numpy array)
#	Output:
#		basis = Value of basis function
#	"""
#	r = norm(x1 - x2);
#	e = 1200;
#	basis = exp(-(e*r)**2);
#	return basis;

#def RBF_Matrix(nodes):										# Old Function, Took too long
#	""" Function to Generate Matrix for RBF interpolation
#	Input Arguments:
#		nodes = N by 2 array containing interpolation nodes, in the form [x, y] in each row (numpy array)
#	Output:
#		A = Matrix A
#	"""
#	print("Generating Matrix for RBF interpolation");
#	ts = time.time();
#	N = len(nodes);
#	A = zeros((N, N));
#	for i in range(N):
#		for j in range(N):
#			A[i, j] = RBF(nodes[i, :], nodes[j, :]);
#	print("That took: ", time.time() - ts);
#	return A;

#def Inter(a, x, y, nodes):									# Old Function, Took too long
#	""" Function to evaluate RBF at desired location
#	Input Arguments:
#		a = Coefficient of each basis function (numpy array)
#		x = x coordinate (numpy array)
#		y = y coordinate (numpy array)
#		nodes = Interpolation nodes
#	Output:
#		fi = Value of RBF at all x, y 's
#	"""
#	print("Reconstructing Function from RBF");
#	ts = time.time();
#	fi = zeros((len(x), len(y)));
#	for i in range(len(x)):
#		for j in range(len(y)):
#			f1 = 0;
#			for g in range(len(nodes)):
#				f1 += a[g]*RBF(array([x[i], y[j]]), nodes[g, :]);
#			fi[i, j] = f1;
#	print("That took: ", time.time() - ts);
#	return fi.T;

def RBF(r):
	""" Function to define Radial Basis Function
	Input Arguments:
		r = distance between two nodes (float or numpy array)
	Output:
		basis = Value of basis function
	"""
	e = 1200;
	basis = exp(-(e*r)**2);
	return basis;

def RBF_Matrix(nodes):
	""" Function to Generate Matrix for RBF interpolation
	Input Arguments:
		nodes = N by 2 array containing interpolation nodes, in the form [x, y] in each row (numpy array)
	Output:
		A = Matrix A
	"""
	print("Generating Matrix for RBF interpolation");
	ts = time.time();
	Norm = lambda x, y: sqrt(x**2 + y**2);
	X, Y = meshgrid(nodes[:, 0], nodes[:, 1]);
	dX = X.T - X;
	dY = Y - Y.T;
	A = RBF(Norm(dX, dY));
	print("That took: ", time.time() - ts);
	return A;

def Inter(a, x, y, nodes):
	""" Function to evaluate RBF at desired location
	Input Arguments:
		a = Coefficient of each basis function (numpy array)
		x = x coordinate (numpy array)
		y = y coordinate (numpy array)
		nodes = Interpolation nodes
	Output:
		fi = Value of RBF at all x, y 's
	"""
	print("Reconstructing Function from RBF");
	ts = time.time();
	Norm = lambda x, y: sqrt(x**2 + y**2);
	ts = time.time();
	fi = zeros((len(x), len(y)));
	for i in range(len(x)):
		for j in range(len(y)):
			r = norm(array([x[i], y[j]]) - nodes, axis = 1);
			fi[i, j] = a.dot(RBF(r));
	print("That took: ", time.time() - ts);
	return fi.T;

def f(x, y): return (1 - x**2)*(1 - y**2);

#%%
# Input Data
chord = 0.515;
span = 2.691;
Nx = 41;
Nz = 81;
x = Coordinate(Nx, span, "x");
z = Coordinate(Nz, chord, "z");
X, Z = meshgrid(x, z);
divx = 1;
divz = 1;
X1, Z1 = meshgrid(x[::divx], z[::divz]);
data = genfromtxt("aerodynamicloaddo228.dat", delimiter = ',');


fig = plt.figure();
ax = plt.axes(projection = '3d');
ax.plot_surface(X, Z, data, rstride = 1, cstride = 1,
				cmap = 'jet', edgecolor = 'none');
#plt.plot(X, Z, 'x');
ax.set_xlabel('Span [m]');
ax.set_ylabel('Chord [m]');
ax.set_zlabel('Loading [kN/m^2]');

#%% Interpolation
Nx1 = len(X1[0, :]);
Nz1 = len(Z1[:, 0]);
nodes = array([X1.reshape(1, -1)[0], Z1.reshape(1, -1)[0]]).T;
A = RBF_Matrix(nodes);
if det(A) == 0.0:
	print("Ill Conditioned Problem!");
else:
	b = data[::divx, ::divz].reshape(1, -1)[0];
	coeff = inv(A).dot(b);
	F = Inter(coeff, x, z, nodes);
	#%% Plots
	fig = plt.figure();
	ax = plt.axes(projection = '3d');
	ax.plot_surface(X, Z, F, rstride = 1, cstride = 1,
					cmap = 'jet', edgecolor = 'none');
	plt.plot(X1, Z1, 'x');
	ax.set_xlabel('Span [m]');
	ax.set_ylabel('Chord [m]');
	ax.set_zlabel('Loading [kN/m^2]');
	plt.show();

	#%% Quadrature Integral
	Dop = 11;
	x_int = linspace(0, span, Dop);
	z_int = linspace(-chord, 0, Dop);
	wx = Quadrature_weights(x_int);
	wz = Quadrature_weights(z_int);
	force_test = integrate(ones((Dop, Dop)), wx, wz);
	force = integrate(Inter(coeff, x_int, z_int, nodes), wx, wz);
	print("Aerodynamic force =", force, "[kN]");
	print("test case:");
	print("Analytical Test =", chord*span);
	print("Quadrature =", force_test);