from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import inv, solve, eigvals, det, norm, cond
import QuadratureIntegrals as quadInt
import AeroLoadInterpolation as AeroInterpol

#%% Input Data
chord = 0.515;
span = 2.691;
Nx = 41;
Nz = 81;
x = AeroInterpol.Coordinate(Nx, span, "x");
z = AeroInterpol.Coordinate(Nz, chord, "z");
X, Z = meshgrid(x, z);
data = genfromtxt("aerodynamicloaddo228.dat", delimiter = ',');

#%% Plot Raw Data
fig = plt.figure();
ax = plt.axes(projection = '3d');
ax.plot_surface(X, Z, data, rstride = 1, cstride = 1,
				cmap = 'jet', edgecolor = 'none');
#plt.plot(X, Z, 'x');
ax.set_xlabel('Span [m]');
ax.set_ylabel('Chord [m]');
ax.set_zlabel('Loading [kN/m^2]');

#%% Interpolation
nodes = array([X.reshape(1, -1)[0], Z.reshape(1, -1)[0]]).T;
A = AeroInterpol.RBF_Matrix(nodes);
if det(A) == 0.0:
	print("Ill Conditioned Problem!");
else:
	b = data.reshape(1, -1)[0];
	coeff = inv(A).dot(b);
	x_p = quadInt.Cosine_Sampler(x[0], x[-1], Nx);
	z_p = linspace(z[0], z[-1], Nz);
	X_p, Z_p = meshgrid(x_p, z_p);
	F = AeroInterpol.Inter(coeff, x_p, z_p, nodes);
	err = AeroInterpol.InerpolationError(coeff, A);
	#%% Plots
	fig = plt.figure();
	ax = plt.axes(projection = '3d');
	ax.plot_surface(X_p, Z_p, F, rstride = 1, cstride = 1,
					cmap = 'jet', edgecolor = 'none');
	#plt.plot(X1, Z1, 'x');
	ax.set_xlabel('Span [m]');
	ax.set_ylabel('Chord [m]');
	ax.set_zlabel('Loading [kN/m^2]');

#%% Quadrature Integral
Dop = 16;
x_int_std = linspace(-1, 1, Dop);
w_std = quadInt.Quadrature_weights(x_int_std);
x_int = quadInt.Cosine_Sampler(0, span, Dop);
z_int = quadInt.Cosine_Sampler(-chord, 0, Dop);
wx = quadInt.Quadrature_weights(x_int);
wz = quadInt.Quadrature_weights(z_int);
force_test = quadInt.integrate(ones((Dop, Dop)), stack((wx, wz)), "2D");
force = quadInt.integrate(AeroInterpol.Inter(coeff, x_int, z_int, nodes), stack((wx, wz)), "2D");
f = lambda x, y: AeroInterpol.Inter(coeff, x, y, nodes);
moment = quadInt.Triple_indefIntegral([0, span], [-chord, 0], Dop, f, x_int_std, w_std);
print("Aerodynamic force =", force, "[kN]");
print("Aerodynamic moment =", moment, "[kNm]");
plt.show();