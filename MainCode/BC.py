from numpy import*
#%% Input data
x1 = 0.174;
x2 = 1.051;
x3 = 2.512;
xa = 0.3;
xa1 = x2 - xa/2;
xa2 = x2 + xa/2;

def step(x, p):
	if x <= 0:
		return 0;
	else:
		return x**p;

def theta_twist_bc(x, integral_val, G, J, z_sc, theta_A, p, R, BC):
	theta_bc = zeros((1, 12))[0];
	theta_bc[0] = z_sc*step(x - x1, 1);
	theta_bc[3] = z_sc*step(x - x2, 1);
	theta_bc[5] = z_sc*step(x - x3, 1);
	theta_bc[2] = -sin(theta_A)*(R - z_sc)*step(x - xa1, 1) + cos(theta_A)*R*step(x - xa1, 1);
	theta_bc[11] = 1;
	rhs = -p*sin(theta_A)*(R - z_sc)*step(x - xa2, 1)/(G*J) + p*cos(theta_A)*R*step(x - xa2, 1)/(G*J) - integral_val/(G*J);
	return theta_bc/(G*J), rhs + BC;

def Torque(x, integral_val, z_sc, theta_A, p, R, BC):
	theta_bc = zeros((1, 12))[0];
	theta_bc[0] = z_sc*step(x - x1, 0);
	theta_bc[3] = z_sc*step(x - x2, 0);
	theta_bc[5] = z_sc*step(x - x3, 0);
	theta_bc[2] = -sin(theta_A)*(R - z_sc)*step(x - xa1, 0) + cos(theta_A)*R*step(x - xa1, 0);

	rhs = -p*sin(theta_A)*(R - z_sc)*step(x - xa2, 0) + p*cos(theta_A)*R*step(x - xa2, 0) - integral_val;
	return theta_bc, rhs + BC;

def v_bc(x, p, Izz, E, theta_A, integral_val, BC):
	v = zeros((1, 12))[0];
	v[0] = step(x - x1, 3);
	v[2] = step(x - xa1, 3)*sin(theta_A);
	v[3] = step(x - x2, 3);
	v[5] = step(x - x3, 3); 
	v[7] = 6*x;
	v[8] = 6;

	rhs = p*sin(theta_A)*step(x - xa2, 3) + 6*integral_val;
	
	return v/(6*E*Izz), rhs/(6*E*Izz) + BC;


def w_bc(x, p, Iyy, E, theta_A, integral_val, BC):
	w = zeros((1, 12))[0];
	w[1] = step(x - x1, 3);
	w[2] = cos(theta_A)*step(x - xa1, 3);
	w[4] = step(x - x2, 3);
	w[6] = step(x - x3, 3);
	w[9] = 6*x;
	w[10] = 6;

	rhs = p*cos(theta_A)*step(x - xa2, 3);
	
	return w/(6*Iyy*E), rhs/(6*Iyy*E) + BC;
	

	
def moment_z_shear_y (x, p, theta_A, integral_val, bool, BC):
	if bool == "Moment":
		i = 1;
		rhs = -p*sin(theta_A)*step(x - xa2, i) - integral_val;
		
	elif bool == "Shear":
		i = 0;
		rhs = -p*sin(theta_A)*step(x - xa2, i) - integral_val;
		
	v = zeros((1, 12))[0];
	v[0] = step(x - x1, i);
	v[2] = step(x - xa1, i)*sin(theta_A);
	v[3] = step(x - x2, i);
	v[5] = step(x - x3, i);

	return v, -rhs + BC;

def moment_y_shear_z (x, p, theta_A, integral_val, bool, BC):
	if bool == "Moment":
		i = 1;
		rhs = -p*cos(theta_A)*step(x - xa2, i);
		
	elif bool == "Shear":
		i = 0;
		rhs = -p*cos(theta_A)*step(x - xa2, i);
		
	w = zeros((1, 12))[0];
	w[1] = step(x - x1, i);
	w[2] = cos(theta_A)*step(x - xa1, i);
	w[4] = step(x - x2, i);
	w[6] = step(x - x3, i);

	return w, -rhs + BC;

def deflection_v(x, R1y, Ra1, theta_a, R2y, R3y, P, integral_val, C1, C2, E, Izz, bool):
	if bool == "Moment":
		c = -1;
		i = 1;
		C1 = 0; C2 = 0;
	elif bool == "DEF":
		c = 6*E*Izz;
		integral_val = 6*integral_val;
		i = 3;
	elif bool == "Shear":
		c = -1;
		i = 0;
		C1 = 0; C2 = 0;
	v = 1/c*(R1y*step(x - x1, i) + Ra1*sin(theta_a)*step(x - xa1, i) + R2y*step(x - x2, i) + 
				R3y*step(x - x3, i) - P*sin(theta_a)*step(x - xa2, i) - integral_val + 6*C1*x + 6*C2);
	return v;

def deflection_w(x, R1z, Ra1, theta_a, R2z, R3z, P, integral_val, C1, C2, E, Iyy, bool):
	if bool == "Moment":
		c = -1;
		i = 1;
		C1 = 0; C2 = 0;
	elif bool == "DEF":
		c = 6*E*Iyy;
		integral_val = 6*integral_val;
		i = 3;
	elif bool == "Shear":
		c = -1;
		i = 0;
		C1 = 0; C2 = 0;
	w = 1/c*(R1z*step(x - x1, i) + Ra1*cos(theta_a)*step(x - xa1, i) + R2z*step(x - x2, i) + 
				R3z*step(x - x3, i) - P*cos(theta_a)*step(x - xa2, i) - integral_val + 6*C1*x + 6*C2);
	return w;

def twist(x, R1y, Ra1, theta_a, R2y, R3y, P, integral_val, C, G, J, r, z_sc, bool):
	if bool == "Theta":
		i = 1;
		c = G*J;
	elif bool == "Torque":
		i = 0;
		c = 1;
		C = 0;
	theta = 1/c*(R1y*z_sc*step(x - x1, i) + R2y*z_sc*step(x - x2, i) + R3y*z_sc*step(x - x3, i) - 
			  Ra1*sin(theta_a)*(r - z_sc)*step(x - xa1, i) + Ra1*cos(theta_a)*r*step(x - xa1, i) + 
			  P*sin(theta_a)*(r - z_sc)*step(x - xa2, i) - P*cos(theta_a)*r*step(x - xa2, i) + integral_val + C);
	return theta;