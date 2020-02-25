from numpy import*
from numpy.linalg import inv, det, norm
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
			pos[i - 1] = c*0.5*(l_ref/2*(1 - cos(theta_i)) + l_ref/2*(1 - cos(theta_j)));
	return pos;

## Radial Basis Function Interpolation
def RBF(r):
	""" Function to define Radial Basis Function
	Input Arguments:
		r = distance between two nodes (float or numpy array)
	Output:
		basis = Value of basis function
	"""
	e = 260;
	basis = sqrt(1 + (e*r)**2);
	#basis = exp(-(800*r)**2)
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
	print("That took: ", time.time() - ts, "s");
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
	if shape(x) != (len(x),):
		x = x[0, :];
		y = y[:, 0];
	Norm = lambda x, y: sqrt(x**2 + y**2);
	ts = time.time();
	fi = zeros((len(x), len(y)));
	for i in range(len(x)):
		for j in range(len(y)):
			r = norm(array([x[i], y[j]]) - nodes, axis = 1);
			fi[i, j] = a.dot(RBF(r));
	return fi.T;

def InerpolationError(a, A):
	""" Function to evaluate the interpolation error using the formula
		from litrature : e_i = a_i/inv(A_{i, i})
	Input Arguments:
		a = array of coefficients of each RBF (numpy array or list)
		A = Interpolation matrix (N X N matrix, numpy array)
	Output:
		e = error (numpy array)
	"""
	e = a/diagonal(inv(A));
	print("Interpolation residuals =", e);
	print("Interpolation error =", amax(e));
	return e;