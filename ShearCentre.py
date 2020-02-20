from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm, inv

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
	zero = zeros_like(z)[1: -1];
	y_spar = linspace(-height, height, len(z))[1: -1];
	return y, y1, zero, y_spar;

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

def integrate(f, w):
	""" Function to Compute integral using Quadrature rule
	Input Arguments:
		f = Matrix containing function values at the nodes (numpy array)
		w = Quadrature weights (numpy array)
	Output:
		I = Value of Integral
	"""
	I = f.dot(w);
	return I;

def Quadrature_weightTransform(nodes, weights, a, b):
	""" Function to transform Quadrature weights and nodes from standard domain [-1, 1] to another domain [a, b]
		while maintaining Degree of Percision
	Input Arguments:
		nodes = array with nodes in [-1, 1] (numpy array)
		weights = array with weights calculated for nodes in [-1, 1] (numpy array)
		[a, b] = new interval (float, float)
	Output:
		x_int = new Quadrature nodes (numpy array)
		w = new Quadrature weights (numpy array)
	"""
	x_int = (b - a)/2*(nodes + 1) + a;
	w = (b - a)/2*weights;
	return x_int, w;

def ShearFlow(BC, t, ds, y, s_stringer, q_boom):
	""" Function to compute the shear flow distribution over y, integrated over the differential
		element ds
	Input Arguments:
		BC = Shear Flow boundary condition at a known point in the structure (float) [N/m]
		t = thickness of the structure (float) [m]
		ds = differntial element in the direction of s (numpy array)
	Output:
		qb_array = array containing the value of the shear flows (numpy array) [N/m]
	"""
	qb_array = zeros(len(y));
	qb_array[0] = BC;
	s = zeros_like(y);
	c = 0;
	for i in range(1, len(qb_array)):
		qb_array[i] = qb_array[i - 1] - 1/Izz*t*y[i]*ds[i - 1];
		s[i] = s[i - 1] + ds[i - 1];
		if s[i] >= s_stringer and c == 0:
			qb_array[i] += q_boom;
			print(q_boom)
			c += 1;
	return qb_array, s;

def ShearFlowIntegrate(q, DOP, s):
	""" Function to integrate the shear flow distribution to find shear force
	Input Arguments:
		q = shear flow distribution (numpy array) [N/m]
		DOP = Degree of Precision of the Quadrature rule desired (int)
		s = vector of s coordinates (numpy array)
	Output:
		Sf = Shear Force (float) [N]
	"""
	global x_int_std, w_std;
	idx_skip = int(ceil(len(s)/DOP));
	q = q[::idx_skip]; s = s[::idx_skip];
	x_int, w = Quadrature_weightTransform(x_int_std, w_std, s[0], s[-1]);
	Sf = integrate(q, w);
	return Sf;


#%% Input Data
N = 10;
h = 0.248; r = h/2;
t = 1.1e-3;
t_spar = 2.2e-3;
Izz = 1.42215e-5;
length = 0.515;
theta = linspace(0, pi/2, N);
z_cord = linspace(0, length - r, N);
z_arc1, y_arc1 = Arc(h/2, theta);				# Arc 1 = upper arc
y_arc2 = -y_arc1;								# Arc 2 = lower arc
y_, y_1, zero, y_spar = Triangel(h/2, length - r, z_cord);
DOP = 3;
x_int_std = linspace(-1, 1, DOP);
w_std = Quadrature_weights(x_int_std);
A_stringer = 5.4e-1;
stringer_no = 11;
stringer_pos = array([[0.1215 , 0.0 ], 
[0.07675071331511048 , 0.09418647580945623  ],
[-0.02327089923050003 , 0.11399894755861381 ],
[-0.12811921373607146 , 0.08074781968472415 ],
[-0.2329675282416429 , 0.047496691810834483 ],
[-0.3378158427472144 , 0.014245563936944823 ],
[-0.3378158427472142 , -0.01424556393694485 ],
[-0.23296752824164274 , -0.04749669181083451],
[-0.1281192137360713 , -0.08074781968472418 ],
[-0.023270899230499863 , -0.11399894755861384],
[0.07675071331511048 , -0.09418647580945623 ]]);		# (z, y)

stringer_pos_s = array([0.0, 0.10999458035397217, 
						0.2199891502692524, 0.32998372175408963, 
						0.439978302771062, 0.5499728742558994, 
						0.659967452249893, 0.7699620237347303, 
						0.8799566047517027, 0.9899511762365399, 1.09994574615182]);


#%% Section 1 & 2 (Noa's Arc)
""" This section computes the shear flow in the LE arc of the Aileron
Due to symmetry shear flow in section 1 and 2 are the same but in opposet
direction.
"""
dtheta = absolute(arctan(y_arc1[1:]/z_arc1[1:]) - arctan(y_arc1[: -1]/z_arc1[:-1]));
q_boom = -1/Izz*A_stringer*stringer_pos[1, 1]; 
qb_skin1, s1 = ShearFlow(0, t, r*dtheta, y_arc1, s_[1], q_boom);

qb_skin2 = -qb_skin1;

##%% Section 3 & 4 (Vertical Spar)
#""" This section computes the shear flow in the vertical spar (contributing to section of the arc) of the Aileron
#Due to symmetry shear flow in section 3 and 4 are in the same direction.
#"""
y_spar1 = y_spar[: int(len(y_spar)/2)];			# y coordinates from the bottom of the spar to the mid point
y_spar2 = y_spar[int(len(y_spar)/2) :];			# y coordinates from the mid point to the top

#ds = (y_spar1[1] - y_spar1[0])*ones_like(y_spar1);
#qb_skin3, s3 = ShearFlow(qb_skin2[-1], t_spar, ds, y_spar1);
#qb_skin4 = qb_skin3[::-1];

#%% Section 3 & 4 (Vertical Spar)
""" This section computes the shear flow in the vertical spar (contributing to section of the triangle) of the Aileron
Due to symmetry shear flow in section 1 and 2 are in the same direction. Cut is made at the mid point to preserve symmetry
"""
ds = (y_spar1[1] - y_spar1[0])*ones_like(y_spar1);
qb_skin5, s5 = ShearFlow(0, t_spar, ds, y_spar2, 0, 0);
qb_skin6 = -qb_skin5[::-1];

#%% Section 5 & 6 (Arms of the triangle)
""" This section computes the shear flow in the two arams of the triangle (contributing to section of the triangle) of the Aileron
Due to symmetry shear flow in section 5 and 6 are the same but in opposite directions. Cut is made at the mid point to preserve symmetry
"""
dy = y_[1] - y_[0]; dz = z_cord[1] - z_cord[0];
ds = norm(array([dy, dz]))*ones_like(y_);
qb_skin7, s7 = ShearFlow(qb_skin5[-1] - qb_skin2[-1], t, ds, y_, 0, 0);
qb_skin8 = -qb_skin7;

#%% Verification
Fz = 0;
Fy = 0;
for i in range(len(qb_skin1) - 1):
	theta = arctan(y_arc1[i]/z_arc1[i]);
	dtheta = abs(arctan(y_arc1[i + 1]/z_arc1[i + 1]) - theta);
	ds = r*dtheta;
	Fz -= cos(pi/2 - theta)*qb_skin1[i]*ds;
	Fy += sin(pi/2 - theta)*qb_skin1[i]*ds;


for i in range(len(qb_skin2) - 1):
	theta = arctan(y_arc2[i]/z_arc1[i]);
	dtheta = abs(arctan(y_arc2[i - 1]/z_arc1[i - 1]) - theta);
	ds = r*dtheta;
	Fz -= cos(pi/2 - abs(theta))*qb_skin2[i]*ds;
	Fy -= sin(pi/2 - abs(theta))*qb_skin2[i]*ds;

#for i in range(len(qb_skin3) - 1):
#	ds = y_spar2[i + 1] - y_spar2[i];
#	Fy += 2*qb_skin3[i]*ds;
#Fy += 2*ShearFlowIntegrate(qb_skin3, DOP, s3);
	
for i in range(len(qb_skin5) - 1):
	ds = y_spar2[i + 1] - y_spar2[i];
	Fy += qb_skin5[i]*ds;
	
for i in range(len(qb_skin6) - 1):
	ds = y_spar1[i + 1] - y_spar1[i];
	Fy -= qb_skin6[i]*ds;

for i in range(len(qb_skin7) - 1):
	dy = y_[i + 1] - y_[i]; dz = z_cord[i + 1] - z_cord[i];
	ds = norm(array([dy, dz]));
	theta = abs(arctan(dy/dz));
	Fz -= cos(theta)*qb_skin7[i]*ds;
	Fy -= sin(theta)*qb_skin7[i]*ds;

for i in range(len(qb_skin8) - 1):
	dy = y_[i + 1] - y_[i]; dz = z_cord[i + 1] - z_cord[i];
	ds = norm(array([dy, dz]));
	theta = abs(arctan(dy/dz));
	Fz -= cos(theta)*qb_skin8[i]*ds;
	Fy += sin(theta)*qb_skin8[i]*ds;

print("Fz =", Fz);
print("Fy =", Fy);

#%% Plots
fig = plt.figure(figsize = (8, 8));
ax = plt.axes(projection = '3d');
# Arc Shear Flows
ax.scatter(-z_arc1, y_arc1, qb_skin1);
ax.scatter(-z_arc1, y_arc2, qb_skin2);
# Spar Shear Flows
#ax.scatter(zero[: int(len(zero)/2)], y_spar1, qb_skin3);		# Section (Arc)
ax.scatter(zero[: int(len(zero)/2)], y_spar1, qb_skin6);		# Section (Spar)
#ax.scatter(zero[int(len(zero)/2) :], y_spar2, qb_skin4);		# Section (Arc)
ax.scatter(zero[int(len(zero)/2) :], y_spar2, qb_skin5);		# Section (Spar)
# Trangle arms shear flow
ax.scatter(z_cord, y_, qb_skin7);
ax.scatter(z_cord, y_1, qb_skin8);
plt.plot(-z_arc1, y_arc1, 'x');
plt.plot(-z_arc1, y_arc2, 'x');
plt.plot(z_cord, y_, 'x');
plt.plot(z_cord, y_1, 'x');
plt.plot(zero, y_spar, 'x');
plt.grid(True);
ax.set_xlabel('Span [m]');
ax.set_ylabel('Chord [m]');
ax.set_zlabel('Shear Flow [N/m]');
plt.show();