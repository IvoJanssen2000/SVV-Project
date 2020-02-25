# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:28:12 2020

@author: elraw
"""

from numpy import *

def theta_twist_bc(x, integral_val, G, J,z_sc,theta_A,p,R):
    theta_bc = zeros((1,12))[0]
    theta_bc[0] = z_sc * step(x - x1)
    theta_bc[3] = z_sc * step(x - x2)
    theta_bc[5] = z_sc * step(x - x3)
    theta_bc[2] = -sin(theta_A) *(R-z_sc)*step(x-xa1)+ cos(theta_A)*R*step(x-xa1)
    
    rhs = -p*sin(theta_A)(R-z_sc)*step(x-xa2) + p*cos(theta_A)*R*step(x-xa2) - integral_val
    return theta_bc/(G*J) , rhs
        
def step(x,theta_A):
    if x<0:
        return(0)
    else:
        return(x)
    
    
    
def v_bc (x,p, Izz,E, theta_A, integral_val):
    
    v = zeros((1,12))[0]
    v[0] = step(x-x1)**3
    v[2] = step(x-xa1)**3 * sin(theta_A)
    v[3] = step(x-x2)**3
    v[5] = step(x-x3)**3 
    v[7] = 6*x
    v[8] = 6

    rhs = p*sin(theta_A) *step(x-xa2)**3 + 6*integral_val
    
    return v/(6*E*Izz), rhs/(6*E*Izz)




def w_bc (x,p, Iyy,E, theta_A, integral_val):
    
    w = zeros((1,12))[0]
    w[1] = step(x-x1)**3 
    w[3] = cos(theta_A)*step(x-xa1)**3
    w[4] = step(x-x2)**3 
    w[6] = step(x-x3)**3
    w[9] = 6*x
    w[10] = 6
    rhs = p*cos(theta_A)* step(x-xa2)**3
    
    return w/(6*Iyy*E), rhs/(6*Iyy*E)
    

    
def moment_z_shear_y (x,p, theta_A, integral_val, bool):
    
    if bool == "Moment":
        i = 1;
        rhs = -p*sin(theta_A) *step(x-xa2)**i - integral_val
        
    elif bool == "Shear":
        i = 0;
        rhs = -p*sin(theta_A) *step(x-xa2)**i - integral_val
        
    v = zeros((1,12))[0]
    v[0] = step(x-x1)**i
    v[2] = step(x-xa1)**i * sin(theta_A)
    v[3] = step(x-x2)**i
    v[5] = step(x-x3)**i 
    
    
    return -v, rhs

def moment_y_shear_z (x,p, theta_A, integral_val, bool):
    
    if bool == "Moment":
        i = 1;
        rhs = -p*cos(theta_A)*step(x-xa2)**i
        
    elif bool == "Shear":
        i = 0;
        rhs = -p*cos(theta_A)*step(x-xa2)**i
        
    w = zeros((1,12))[0]
    w[1] = step(x-x1)**i
    w[3] = cos(theta_A)*step(x-xa1)**i
    w[4] = step(x-x2)**i
    w[6] = step(x-x3)**i
   
    
    return  -w , rhs 
   
    