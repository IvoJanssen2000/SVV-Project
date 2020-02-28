from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import inv, solve, eigvals, det, norm, cond
import QuadratureIntegrals as quadInt
import AeroLoadInterpolation as AeroInterpol
import BC
import time

#%% Input Data Do228
chord = 0.515;
span = 2.691;
x1 = 0.174;
x2 = 1.051;
x3 = 2.512;
xa = 0.3;
xa1 = x2 - xa/2;
xa2 = x2 + xa/2;
J = 1.9193312e-05;
G = 28e9;
Izz = 1.4004e-5;
Iyy = 5.37e-5;
theta_aileron = 25*pi/180;
r = 0.124;
E = 73.1e9;
d1 = 1.034e-2;
d3 = 2.066e-2;
P = 20.6e3;
z_sc = -0.0023667;
Nx = 41;
Nz = 81;
x = AeroInterpol.Coordinate(Nx, span, "x");
z = AeroInterpol.Coordinate(Nz, chord, "z");
X, Z = meshgrid(x, z);
data = genfromtxt("aerodynamicloaddo228.dat", delimiter = ',');
z_vect = [-chord, 0];
# B737
#chord = 0.605;
#span = 2.661;
#x1 = 0.172;
#x2 = 1.211;
#x3 = 2.591;
#xa = 35.0/100;
#xa1 = x2 - xa/2;
#xa2 = x2 + xa/2;
#J = 1.9193312e-05;
#G = 28e9;
#E = 73.1e9;
#r = 20.5/100;
#d1 = 1.154/100;
#d3 = 1.840/100;
#theta_aileron = 28*pi/180;
#P = 97.4*1000;
#Izz = 7.287102054648984e-06;
#Iyy = 2.5203170137972818e-05;
#Nx = 50;
#Nz = 2;
#z_sc = r - 0.24718604786615267;
#x = linspace(0, span, Nx);
#z = linspace(0.25*chord, 0.25*chord + 0.001, Nz);
#X, Z = meshgrid(x, z);
#data = 5.54*ones_like(X);
#z_vect = [0.25*chord, 0.25*chord + 0.001];

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
	f = lambda x, y: AeroInterpol.Inter(coeff, x, y, nodes);
	#%% Plots
	fig = plt.figure();
	ax = plt.axes(projection = '3d');
	ax.plot_surface(X_p, Z_p, F, rstride = 1, cstride = 1,
					cmap = 'jet', edgecolor = 'none');
	#plt.plot(X1, Z1, 'x');
	ax.set_xlabel('Span [m]');
	ax.set_ylabel('Chord [m]');
	ax.set_zlabel('Loading [kN/m^2]');
	print("Done Interpolation");
#%% Quadrature Integral
Dop = 10;
x_int_std = linspace(-1, 1, Dop);
w_std = quadInt.Quadrature_weights(x_int_std);
x_int, wx = quadInt.Quadrature_weightTransform(x_int_std, w_std, 0, span);
z_int, wz = quadInt.Quadrature_weightTransform(x_int_std, w_std, z_vect[0], z_vect[1]);
force = quadInt.QuadratureInt(AeroInterpol.Inter(coeff, x_int, z_int, nodes), stack((wx, wz)), "2D");

#%% Assembly of Matrix
print("Assembling Matrix");
N = 12;
A = zeros((N, N));
b = zeros((1, N))[0];
## BC at x1
Intx0 = 1000*quadInt.FiveD_indefIntegral([0, x1], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std);
theta_intx0 = 1000*quadInt.Triple_indefIntegral([0, x1], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, z_sc);
A[0, :], b[0] = BC.v_bc(x1, P, Izz, E, theta_aileron, Intx0, d1*cos(theta_aileron));
A0, b0 = BC.theta_twist_bc(x1, theta_intx0, G, J, z_sc, theta_aileron, P, r, 0);
A[0, :] += A0*z_sc; b[0] += b0*z_sc;
A[1, :], b[1] = BC.w_bc(x1, P, Iyy, E, theta_aileron, 0, -d1*sin(theta_aileron));

## BC at x2
Intx1 = 1000*quadInt.FiveD_indefIntegral([0, x2], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std);
theta_intx2 = 1000*quadInt.Triple_indefIntegral([0, x2], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, z_sc);
A[2, :], b[2] = BC.v_bc(x2, P, Izz, E, theta_aileron, Intx1, 0);
A2, b2 = BC.theta_twist_bc(x2, theta_intx2, G, J, z_sc, theta_aileron, P, r, 0);
A[2, :] += A2*z_sc; b[2] += b2*z_sc;
A[3, :], b[3] = BC.w_bc(x2, P, Iyy, E, theta_aileron, 0, 0);

## BC at x3
Intx3 = 1000*quadInt.FiveD_indefIntegral([0, x3], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std);
theta_intx3 = 1000*quadInt.Triple_indefIntegral([0, x3], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, z_sc);
A[4, :], b[4] = BC.v_bc(x3, P, Izz, E, theta_aileron, Intx3, d3*cos(theta_aileron));
A4, b4 = BC.theta_twist_bc(x3, theta_intx3, G, J, z_sc, theta_aileron, P, r, 0);
A[4, :] += A4*z_sc; b[4] += b4*z_sc;
A[5, :], b[5] = BC.w_bc(x3, P, Iyy, E, theta_aileron, 0, -d3*sin(theta_aileron));

## BC at x = la
W = stack((wx, wz));
XX, ZZ = meshgrid(x_int, z_int);
F = f(XX, ZZ);
A[6, :], b[6] = BC.moment_y_shear_z(span, P, theta_aileron, 0, "Moment", 0);
A[7, :], b[7] = BC.moment_z_shear_y(span, P, theta_aileron,
				1000*quadInt.Triple_indefIntegral([0, span], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, 0), "Moment", 0);
A[8, :], b[8] = BC.moment_y_shear_z(span, P, theta_aileron, 0, "Shear", 0);
A[9, :], b[9] = BC.moment_z_shear_y(span, P, theta_aileron, 1000*quadInt.QuadratureInt(F, W, "2D"), "Shear", 0);
A[10, :], b[10] = BC.Torque(span, 1000*quadInt.QuadratureInt(F*(z_int - z_sc).reshape(len(z_int), 1), W, "2D"), z_sc, theta_aileron, P, r, 0);

## BC at xa1
Intxa1 = 1000*quadInt.FiveD_indefIntegral([0, xa1], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std);
theta_intxa1 = 1000*quadInt.Triple_indefIntegral([0, xa1], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, z_sc);
va1, ba1 = BC.v_bc(xa1, P, Izz, E, theta_aileron, Intxa1, 0);
theta_xa1, ba2 = BC.theta_twist_bc(xa1, theta_intxa1, G, J, z_sc, theta_aileron, P, r, 0);
wa1, ba3 = BC.w_bc(xa1, P, Iyy, E, theta_aileron, 0, 0);
A[11, :] = (va1 - theta_xa1*(r - z_sc))*sin(theta_aileron) + (wa1 + theta_xa1*r)*cos(theta_aileron);
b[11] = ba1*sin(theta_aileron) + ba2*(-(r - z_sc)*sin(theta_aileron) + r*cos(theta_aileron)) + ba3*cos(theta_aileron);
print("Matrix Assembled, Solving for Unknows");
#%%
Soln = inv(A).dot(b);
Ry_total = Soln[0] + Soln[2]*sin(theta_aileron) + Soln[3] + Soln[5];
Fy = P*sin(theta_aileron) + 1000*force;
print("reaction force in y =", Ry_total);
print("Fy =", Fy);
Rz_total = Soln[1] + Soln[2]*cos(theta_aileron) + Soln[4] + Soln[6];
Fz = P*cos(theta_aileron);
print("reaction force in z =", Rz_total);
print("Fz =", Fz);

#%% Deflection
print("Calculating Deflections and Loading Distributions");
ts = time.time();
x = linspace(0, span, 20);
#x = array([2591.       , 1386.       , 2661.       , 1211.       ,
#	   1036.       ,  172.       ,    0.       , 2565.89575  ,
#	   2540.79175  , 2515.6875   , 2490.58325  , 2465.47925  ,
#	   2440.375    , 2415.27075  , 2390.16675  , 2365.0625   ,
#	   2339.95825  , 2314.85425  , 2289.75     , 2264.64575  ,
#	   2239.54175  , 2214.4375   , 2189.33325  , 2164.22925  ,
#	   2139.125    , 2114.02075  , 2088.91675  , 2063.8125   ,
#	   2038.70837  , 2013.60413  , 1988.5      , 1963.39587  ,
#	   1938.29163  , 1913.1875   , 1888.08337  , 1862.97913  ,
#	   1837.875    , 1812.77087  , 1787.66663  , 1762.5625   ,
#	   1737.45837  , 1712.35413  , 1687.25     , 1662.14587  ,
#	   1637.04163  , 1611.9375   , 1586.83337  , 1561.72913  ,
#	   1536.625    , 1511.52087  , 1486.41663  , 1461.3125   ,
#	   1436.20837  , 1411.10413  , 2614.33325  , 2637.66675  ,
#	   1186.       , 1161.       , 1136.       , 1111.       ,
#	   1086.       , 1061.       , 1361.       , 1336.       ,
#	   1311.       , 1286.       , 1261.       , 1236.       ,
#		147.428574 ,  122.85714  ,   98.2857132,   73.7142868,
#		 49.1428566,   24.5714283, 1011.31427  ,  986.628601 ,
#		961.942871 ,  937.257141 ,  912.571411 ,  887.885742 ,
#		863.200012 ,  838.514282 ,  813.828552 ,  789.142883 ,
#		764.457153 ,  739.771423 ,  715.085693 ,  690.400024 ,
#		665.714294 ,  641.028564 ,  616.342834 ,  591.657166 ,
#		566.971436 ,  542.285706 ,  517.599976 ,  492.914276 ,
#		468.228577 ,  443.542847 ,  418.857147 ,  394.171417 ,
#		369.485718 ,  344.799988 ,  320.114288 ,  295.428558 ,
#		270.742859 ,  246.057144 ,  221.371429 ,  196.685715 ]);
#x /= 1000;
Mz = zeros_like(x);
w = zeros_like(x);
Sz = zeros_like(x);
My = zeros_like(x);
v = zeros_like(x);
Sy = zeros_like(x);
T = zeros_like(x);
twist = zeros_like(x);
for i in range(len(x)):
	My[i] = BC.deflection_w(x[i], Soln[1], Soln[2], theta_aileron, Soln[4], Soln[6], P, 0, Soln[9], Soln[10], E, Iyy, "Moment");
	w[i] = BC.deflection_w(x[i], Soln[1], Soln[2], theta_aileron, Soln[4], Soln[6], P, 0, Soln[9], Soln[10], E, Iyy, "DEF");
	Sz[i] = BC.deflection_w(x[i], Soln[1], Soln[2], theta_aileron, Soln[4], Soln[6], P, 0, Soln[9], Soln[10], E, Iyy, "Shear");

	Mz[i] = BC.deflection_v(x[i], Soln[0], Soln[2], theta_aileron, Soln[3], Soln[5], P, 
						 1000*quadInt.Triple_indefIntegral([0, x[i]], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, 0), Soln[7], Soln[8], E, Izz, "Moment");
	v[i] = BC.deflection_v(x[i], Soln[0], Soln[2], theta_aileron, Soln[3], Soln[5], P, 
						 1000*quadInt.FiveD_indefIntegral([0, x[i]], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std), Soln[7], Soln[8], E, Izz, "DEF");
	x_int, w_x = quadInt.Quadrature_weightTransform(x_int_std, w_std, 0, x[i]);
	W = stack((w_x, wz));
	XX, ZZ = meshgrid(x_int, z_int);
	F = f(XX, ZZ);
	Sy[i] = BC.deflection_v(x[i], Soln[0], Soln[2], theta_aileron, Soln[3], Soln[5], P, 
						 1000*quadInt.QuadratureInt(F, W, "2D"), Soln[7], Soln[8], E, Izz, "Shear");
	T[i] = BC.twist(x[i], Soln[0], Soln[2], theta_aileron, Soln[3], Soln[5], P, 
				 1000*quadInt.QuadratureInt(F*(z_int - z_sc).reshape(len(z_int), 1), W, "2D"), Soln[11], G, J, r, z_sc, "Torque");
	twist[i] = BC.twist(x[i], Soln[0], Soln[2], theta_aileron, Soln[3], Soln[5], P, 
				 1000*quadInt.Triple_indefIntegral([0, x[i]], [z_vect[0], z_vect[1]], Dop, f, x_int_std, w_std, z_sc), Soln[11], G, J, r, z_sc, "Theta");
print("That took: ", time.time() - ts, "s");
#%%
fig, ax = plt.subplots(1, 3, figsize = (12, 4));
ax[0].plot(x, My);
ax[0].set_title('Moment in y direction');
ax[0].grid(True);
ax[0].set_xlabel("x [m]");
ax[0].set_ylabel("My [Nm]");
plt.tight_layout();

ax[1].plot(x, w);
ax[1].set_title('Deflection in z direction');
ax[1].grid(True);
ax[1].set_xlabel("x [m]");
ax[1].set_ylabel("w [m]");
plt.tight_layout();

ax[2].plot(x, Sz);
ax[2].set_title('Shear force in z direction');
ax[2].grid(True);
ax[2].set_xlabel("x [m]");
ax[2].set_ylabel("Sz [N]");
plt.tight_layout();

fig, ax = plt.subplots(1, 3, figsize = (12, 4));
ax[0].plot(x, Mz);
ax[0].set_title('Moment in z direction');
ax[0].grid(True);
ax[0].set_xlabel("x [m]");
ax[0].set_ylabel("Mz [Nm]");
plt.tight_layout();

ax[1].plot(x, v + twist*z_sc);
ax[1].set_title('Deflection in y direction');
ax[1].grid(True);
ax[1].set_xlabel("x [m]");
ax[1].set_ylabel("v [m]");
plt.tight_layout();

ax[2].plot(x, Sy);
ax[2].set_title('Shear force in y direction');
ax[2].grid(True);
ax[2].set_xlabel("x [m]");
ax[2].set_ylabel("Sy [N]");
plt.tight_layout();

fig, ax = plt.subplots(1, 2, figsize = (8, 4));
ax[0].plot(x, -T);
ax[0].grid(True);
ax[0].set_title("Torque distribution");
ax[0].set_xlabel("x [m]");
ax[0].set_ylabel("T [Nm]");
plt.tight_layout();

ax[1].plot(x, -twist);
ax[1].grid(True);
ax[1].set_title("Twist");
ax[1].set_xlabel("x [m]");
ax[1].set_ylabel("θ [rad]");
plt.tight_layout();
plt.show();
