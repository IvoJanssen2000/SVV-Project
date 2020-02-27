from numpy import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm, inv, det

def Arc(r, theta):
	""" Function to compute the y and z coordinates of the LE arc 
	Input Arguments:
		r = radius of the arc (float) [m]
		theta = range of theta values (numpy array) [rad]
	Output:
		y, z coordinates (numpy array)
	"""
	y = r*sin(theta);
	z = r*cos(theta);
	return z, y;

def Triangel(height, chord, z):
	""" Function to compute the y and z coordinates of the trangel shape 
	Input Arguments:
		h = Height of the triangle (float) [m]
		chord = chord of the airfoil - radius of arc (float) [m]
	Output:
		y, z coordinates (numpy array)
	"""
	dydz = -height/chord;
	y = height + dydz*z;
	y1 = -height - dydz*z;
	zero = zeros_like(z);
	y_spar = linspace(-height, height, len(z));
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

def ShearFlow(BC, t, Izz, ds, y, s_stringer, q_boom,S):
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
		qb_array[i] = qb_array[i - 1] - S/Izz*t*y[i]*ds[i - 1];
		s[i] = s[i - 1] + ds[i - 1];
		check = len(s_stringer)*[True];
		if (c < len(s_stringer)):
			if s[i] >= s_stringer[c] and check[c]:
				qb_array[i] += q_boom[c];
				check[c] = False;
				c += 1;
	return qb_array, s;

def ShearFlowIntegrate(q, ds, t):
	""" Function to integrate the shear flow distribution to find shear force
	Input Arguments:
		q = shear flow distribution (numpy array) [N/m]
		ds = vector of steps in s(numpy array)
	Output:
		Sf = Shear Force (float) [N]
	"""
	Sf = q.dot(ds)/t;
	return Sf;

def Verification(qb_skin):
	Fz = 0;
	Fy = 0;
	for i in range(len(qb_skin[0]) - 1):
		theta = arctan(y_arc1[i]/z_arc1[i]);
		dtheta = abs(arctan(y_arc1[i + 1]/z_arc1[i + 1]) - theta);
		ds = r*dtheta;
		Fz -= cos(pi/2 - theta)*qb_skin[0][i]*ds;
		Fy += sin(pi/2 - theta)*qb_skin[0][i]*ds;

	for i in range(len(qb_skin[1]) - 1):
		theta = arctan(y_arc2[i]/z_arc1[i]);
		dtheta = abs(arctan(y_arc2[i + 1]/z_arc1[i + 1]) - theta);
		ds = r*dtheta;
		Fz -= cos(pi/2 - abs(theta))*qb_skin[1][i]*ds;
		Fy -= sin(pi/2 - abs(theta))*qb_skin[1][i]*ds;
	
	for i in range(len(qb_skin[2]) - 1):
		ds = y_spar[i + 1] - y_spar[i];
		Fy += qb_skin[2][i]*ds;

	for i in range(len(qb_skin[3]) - 1):
		dy = y_[i + 1] - y_[i]; dz = z_cord[i + 1] - z_cord[i];
		ds = norm(array([dy, dz]));
		theta = abs(arctan(dy/dz));
		Fz -= cos(theta)*qb_skin[3][i]*ds;
		Fy -= sin(theta)*qb_skin[3][i]*ds;

	for i in range(len(qb_skin[4]) - 1):
		dy = y_[i + 1] - y_[i]; dz = z_cord[i + 1] - z_cord[i];
		ds = norm(array([dy, dz]));
		theta = abs(arctan(dy/dz));
		Fz -= cos(theta)*qb_skin[4][i]*ds;
		Fy += sin(theta)*qb_skin[4][i]*ds;
	print("Fz =", Fz, "N");
	print("Fy =", Fy, "N,", "percentage error =", (1 - Fy)*100, "%");




def calculateShear(Sz,Sy,Iyy,Izz):
	#%% Section 1 & 2 (Noa's Arc)
	""" This section computes the shear flow in the LE arc of the Aileron
	Due to symmetry shear flow in section 1 and 2 are the same but in opposet
	direction.
	"""
	dtheta1 = absolute(arctan(y_arc1[1]/z_arc1[1]) - arctan(y_arc1[0]/z_arc1[0]))*ones_like(y_arc1);
	q_boom = -1/Izz*A_stringer*stringer_pos[:, 1]; 
	qb_skin1, s1 = ShearFlow(0, t, Izz, r*dtheta1, y_arc1, [stringer_pos_s[1]], Sy*[q_boom[1]],Sy);
	qb_skin2, _ =  ShearFlow(0, t, Izz, r*dtheta1, y_arc2, [stringer_pos_s[1]], Sy*[-q_boom[1]],Sy); 

	qb_skin1 += ShearFlow(0,t,Iyy,r*dtheta1,z_arc1,[stringer_pos_s[1]],[0],Sz)[0]
	qb_skin2 += ShearFlow(0,t,Iyy,r*dtheta1,z_arc1,[stringer_pos_s[1]],[0],Sz)[0]

	#qb_skin2 = -qb_skin1;

	#%% Section 3 & 4 (Vertical Spar)
	""" This section computes the shear flow in the vertical spar (contributing to section of the triangle) of the Aileron
	Due to symmetry shear flow in section 1 and 2 are in the same direction. Cut is made at the mid point to preserve symmetry
	"""
	ds3 = (y_spar1[1] - y_spar1[0])*ones_like(y_spar);
	qb_skin3, s5 = ShearFlow(0, t_spar, Izz, ds3[: int(len(y_spar1))], y_spar2, [0], [0],Sy);
	#qb_skin4 = -qb_skin3[::-1];
	qb_skin4, _ = ShearFlow(0, t_spar, Izz, ds3[: int(len(y_spar1))], y_spar1[::-1], [0], [0],Sy);

	#%% Section 5 & 6 (Arms of the triangle)
	""" This section computes the shear flow in the two arams of the triangle (contributing to section of the triangle) of the Aileron
	Due to symmetry shear flow in section 5 and 6 are the same but in opposite directions. Cut is made at the mid point to preserve symmetry
	"""
	dy = y_[1] - y_[0]; dz = z_cord[1] - z_cord[0];
	ds5 = norm(array([dy, dz]))*ones_like(y_);
	qb_skin5, s7 = ShearFlow(qb_skin3[-1] - qb_skin2[-1], t, Izz, ds5, y_, stringer_pos_s[2: 6] - pi*r/2, Sy*q_boom[2: 6],Sy);
	#qb_skin6 = -qb_skin5;
	qb_skin6, _ = ShearFlow(-(qb_skin3[-1] - qb_skin2[-1]), t, Izz, ds5, y_1, s_max - stringer_pos_s[6: 10][::-1] - pi/2*r, Sy*q_boom[6: 10][::-1],Sy);

	qb_skin5 += ShearFlow(0,t,Iyy,ds5,-z_cord,stringer_pos_s[2:6] - pi*r/2, [0,0,0,0],Sz)[0]
	qb_skin6 += ShearFlow(0,t,Iyy,ds5,-z_cord,stringer_pos_s[2:6] - pi*r/2, [0,0,0,0],Sz)[0]

	qb_skin34 = concatenate((qb_skin3, -qb_skin4));
	qb_skin = stack((qb_skin1, qb_skin2, qb_skin34, qb_skin5, qb_skin6));
	print("## Verifying Basic Shear Flow Distribution ##");
	Verification(qb_skin);

	#%% qs0 Calculation
	b = zeros(2);
	qb_skin34 = concatenate((-qb_skin3, qb_skin4));
	qbI = stack((qb_skin1, -qb_skin2, qb_skin34));
	qb_skin34_II = concatenate((qb_skin3, -qb_skin4));
	qbII = stack((qb_skin34_II, qb_skin5, -qb_skin6));
	qb = [qbI, qbII];
	for i in range(len(b)):
		for j in range(len(qbI)):
			if i == 0:
				if j == 0 or j == 1:
					ds = r*dtheta1;
					thickness = t;
				else:
					ds = ds3;
					thickness = t_spar;
			elif i == 1:
				if j == 0:
					ds = ds3;
					thickness = t_spar;
				else:
					ds = ds5;
					thickness = t;
			b[i] += ShearFlowIntegrate(qb[i][j, :], ds, thickness);
	b = -b;
	M = array([[h/t_spar + pi*r/t, -h/t_spar], [-h/t_spar, h/t_spar + 2*length_arm/t]]);
	qs0 = inv(M).dot(b);
	qs0[0] = (b[0] - M[0, 1]*qs0[1])/M[0, 0];
	print(qs0);

	#%% Verification qs0
	qb_skin1 = qb_skin1 + qs0[0];
	qb_skin2 = qb_skin2 - qs0[0];
	qb_skin3 = qb_skin3 - qs0[0] + qs0[1];
	qb_skin4 = qb_skin4 + qs0[0] - qs0[1];
	qb_skin5 = qb_skin5 + qs0[1];
	qb_skin6 = qb_skin6 - qs0[1];

	qb_skin34 = concatenate((qb_skin3, -qb_skin4));
	qb_skin = stack((qb_skin1, qb_skin2, qb_skin34, qb_skin5, qb_skin6));
	print("## Verifying Total Shear Flow Distribution ##");
	Verification(qb_skin);
	arm = r*cos(arctan(r/(chord - r)));
	Moment = -r*ShearFlowIntegrate(qb_skin1, r*dtheta1, 1);
	Moment += r*ShearFlowIntegrate(qb_skin2, r*dtheta1, 1);
	Moment -= arm*ShearFlowIntegrate(qb_skin5, ds5, 1);
	Moment += arm*ShearFlowIntegrate(qb_skin6, ds5, 1);
	SC = Moment/Sy;
	print("Shear Centre Location (z, y) =", (SC, 0), "[m]");

	return qb_skin1, qb_skin2, qb_skin3,qb_skin4,qb_skin5,qb_skin6

#%% Input Data
N = 1000;
h = 0.248; r = h/2;
t = 1.1e-3;
t_spar = 2.2e-3;
#Izz = 1.42215e-5;
Iyy = 5.15e-5
Izz = 1.4004e-5;
chord = 0.515;
theta = linspace(0, pi/2, N);
z_cord = linspace(0, chord - r, N);
z_arc1, y_arc1 = Arc(h/2, theta);				# Arc 1 = upper arc
y_arc2 = -y_arc1;								# Arc 2 = lower arc
y_, y_1, zero, y_spar = Triangel(h/2, chord - r, z_cord);
y_spar1 = y_spar[: int(len(y_spar)/2)];			# y coordinates from the bottom of the spar to the mid point
y_spar2 = y_spar[int(len(y_spar)/2) :];			# y coordinates from the mid point to the top
A_stringer = 5.4e-5;
length_arm = sqrt(r**2 + (chord - r)**2);
s_max = pi*r + 2*length_arm;
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

stringer_pos_s = array([0.0, 0.10999390226942096,
					   0.22010769635503435, 0.3300018343791974,
					  0.43998605554511216, 0.5499766891204068,
					 0.6599636373853857, 0.7699542709606804,
					0.8799384921265951, 0.9898326301507581, 1.0999464242363712]);

qb_skin1, qb_skin2, qb_skin3,qb_skin4,qb_skin5,qb_skin6 = calculateShear(0,1000,Iyy,Izz)
#%% Plots
fig = plt.figure(figsize = (8, 8));
ax = plt.axes(projection = '3d');
# Arc Shear Flows
ax.scatter(-z_arc1, y_arc1, qb_skin1);
ax.scatter(-z_arc1, y_arc2, qb_skin2);
# Spar Shear Flows
ax.scatter(zero[: int(len(zero)/2)], y_spar1[::-1], qb_skin4);
ax.scatter(zero[int(len(zero)/2) :], y_spar2, qb_skin3);
# Trangle arms shear flow
ax.scatter(z_cord, y_, qb_skin5);
ax.scatter(z_cord, y_1, qb_skin6);
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
