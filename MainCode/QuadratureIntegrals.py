from numpy import*
from numpy.linalg import inv, norm
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
	N = len(x);
	a = x[0]; b = x[-1];
	A_i = GenMatrix(N, x);
	B = zeros((N, 1));
	for i in range(1, N + 1):
		B[i - 1] = (b**i - a**i)/i;
	w = A_i.dot(B);
	return w.T[0];

def QuadratureInt(f, w, dim):
	""" Function to Compute 1D and 2D integral using Quadrature rule
	Input Arguments:
		f = Matrix containing function values at the nodes (numpy array)
		w = Quadrature weights, if dim = 2D then w has to have 2 rows (numpy array)
		dim = Dimention of definite integral, 1D or 2D (str)
	Output:
		I = Value of Integral
	"""
	if dim == "1D":
		I = f.dot(w);
	elif dim == "2D":
		I = f.dot(w[0, :]).dot(w[1, :]);
	return I;

def Quadrature_weightTransform(nodes, weights, a, b):
	""" Function to transform Quadrature weights and nodes from one domain [c, d] to another [a, b]
		while maintaining Degree of Percision
	Input Arguments:
		nodes = array with nodes in [c, d] (numpy array)
		weights = array with weights calculated for nodes in [c, d] (numpy array)
		[a, b] = new interval (float, float)
	Output:
		x_int = new Quadrature nodes (numpy array)
		w = new Quadrature weights (numpy array)
	"""
	x_int = (b - a)/2*(nodes + 1) + a;
	w = (b - a)/2*weights;
	return x_int, w;

def Triple_indefIntegral(x_vec, z_vec, DOP, f, x_int_std, w_std, torque):
	""" Function to compute triple integral with one indefinite integral and two 
		definite integrals.
	Input Arguments:
		x_vec = array containing the bounds of one of the definite integrals (list or array)
		z_vec = array containing the bounds of the other definite integral (list or array)
	Output:
		I = Value of Integral
	"""
	x_int, w_x = Quadrature_weightTransform(x_int_std, w_std, x_vec[0], x_vec[1]);
	z_int, w_z = Quadrature_weightTransform(x_int_std, w_std, z_vec[0], z_vec[1]);
	g = zeros_like(x_int);
	for i in range(len(g)):
		x_t, w_t = Quadrature_weightTransform(x_int_std, w_std, x_vec[0], x_int[i]);
		X, Z = meshgrid(x_t, z_int);
		fi = f(X, Z);
		if torque != 0:
			fi *= (z_int - torque).reshape(len(z_int), 1);
		w = stack((w_t, w_z));
		g[i] = QuadratureInt(fi, w, "2D");

	I = QuadratureInt(g, w_x, "1D");
	return I;

def FiveD_indefIntegral(x_vec, z_vec, DOP, f, x_int_std, w_std):
	x_int, w_x = Quadrature_weightTransform(x_int_std, w_std, x_vec[0], x_vec[1]);
	G = zeros_like(x_int);
	for i in range(len(x_int)):
		ksi_int, w_ksi = Quadrature_weightTransform(x_int_std, w_std, x_int[0], x_int[i]);
		g = zeros_like(ksi_int);
		for j in range(len(ksi_int)):
			g[j] = Triple_indefIntegral([x_int[0], ksi_int[j]], z_vec, DOP, f, x_int_std, w_std, 0);
		G[i] = QuadratureInt(g, w_ksi, "1D");
	I5D = QuadratureInt(G, w_x, "1D");
	return I5D;

""" Verification """
#def f(x, z): return x + z + x**2 + z**2;

##%% Input
#x0 = 0; x1 = 1;
#z0 = 0; z1 = 1;
#DOP = 4;
#x_int_std = linspace(-1, 1, DOP);
#w_std = Quadrature_weights(x_int_std);
#x_int, w_x = Quadrature_weightTransform(x_int_std, w_std, x0, x1);
#z_int, w_z = Quadrature_weightTransform(x_int_std, w_std, z0, z1);

##%% 2D Integral Verification
#print("## Verification of 2D Integral ## \n");
#X, Z = meshgrid(x_int, z_int);
#fi = f(X, Z);
#w = stack((w_x, w_z));
#I = QuadratureInt(fi, w, "2D");
#analytical = 5/3;
#print("Analytical Solution =", analytical);
#print("Numerical Solution =", I);
#error = (I - analytical)/analytical*100;
#print("Percentage error =", error, "%", "\n");

##%% 3D Integral with indefinite integral Verification
#print("## Verification of 3D Integral with indefinite integral ## \n");
#I = Triple_indefIntegral([x0, x1], [z0, z1], 4, f, x_int_std, w_std, 0);
#analytical = 1/6 + 1/6 + 1/4 + 1/12;
#print("Analytical Solution =", analytical);
#print("Numerical Solution =", I);
#error = (I - analytical)/analytical*100;
#print("Percentage error =", error, "% \n");

##%% 5D Integral with indefinite integral Verification
#print("## Verification of 5D Integral with indefinite integral ## \n");
#G = zeros_like(x_int);
#for i in range(len(x_int)):
#	ksi_int, w_ksi = Quadrature_weightTransform(x_int_std, w_std, x0, x_int[i]);
#	g = zeros_like(ksi_int);
#	for j in range(len(ksi_int)):
#		g[j] = Triple_indefIntegral([x0, ksi_int[j]], [z0, z1], DOP, f, x_int_std, w_std, 0);
#	G[i] = QuadratureInt(g, w_ksi, "1D");

#I5D = QuadratureInt(G, w_x, "1D");
#print(FiveD_indefIntegral([x0, x1], [z0, z1], DOP, f, x_int_std, w_std));
#analytical = 1/120 + 1/48 + 1/72 + 1/360;
#print("Analytical Solution =", analytical);
#print("Numerical Solution =", I5D);
#error = (I5D - analytical)/analytical*100;
#print("Percentage error =", error, "% \n");

##%% RBF integral test
#def RBF(x, y): return sqrt(1 + 400**2*(x**2 + y**2));
#DOP = 15;
#x_int = linspace(-1, 1, DOP);
#z_int = linspace(-1, 1, DOP);
#w_x = Quadrature_weights(x_int);
#w_z = Quadrature_weights(z_int);
#X, Z = meshgrid(x_int, z_int);
#fi = RBF(X, Z);
#fig = plt.figure();
#ax = plt.axes(projection = '3d');
#ax.plot_surface(X, Z, fi, rstride = 1, cstride = 1,
#				cmap = 'jet', edgecolor = 'none');
#w = stack((w_x, w_z));
#I = QuadratureInt(fi, w, "2D");
#analytical = 1224.32;
#print("Analytical Solution =", analytical);
#print("Numerical Solution =", I);
#error = (I - analytical)/analytical*100;
#print("Percentage error =", error, "% \n");
##plt.show();

#I = Triple_indefIntegral([0, 1], [0, 10], DOP, RBF, x_int_std, w_std, 0);
#print(I);