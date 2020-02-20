# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:46:19 2020

@author: marti
"""

import math as m
import numpy as np

Nstringer = 11
Nspar = 1
Npanel = 2
Nsemi = 1
Ca = 0.515      #m
h = 0.248       #m 
tsk = 0.0011    #m
tsp = 0.0022    #m
tst = 0.0012    #m
hst = 0.015     #m
wst = 0.03      #m

R = h/2
r = h/2 - tsk
theta = m.atan(R/(Ca-R))
l = np.sqrt((Ca-R)**2 + R**2)

L = np.pi*R + 2* l
d = L/11
circ = 2*np.pi*R
alpha = m.radians(360 * d/circ)
z = R*m.cos(alpha)
y = R*m.sin(alpha)


def stringerlstgen():
    z6 = -Ca + R + d/2 * m.cos(theta)
    z5 = -Ca + R + d/2 * m.cos(theta) + 1*d*m.cos(theta)
    z4 = -Ca + R + d/2 * m.cos(theta) + 2*d*m.cos(theta)
    z3 = -Ca + R + d/2 * m.cos(theta) + 3*d*m.cos(theta)
    y6 = d/2 * m.sin(theta)
    y5 = d/2 * m.sin(theta) + 1*d*m.sin(theta)
    y4 = d/2 * m.sin(theta) + 2*d*m.sin(theta)
    y3 = d/2 * m.sin(theta) + 3*d*m.sin(theta)
    stringerlst = np.zeros((11,2)) #z,y
    stringerlst[0] = (R,0)
    stringerlst[1] = (z,y)
    stringerlst[2] = (z3,y3)
    stringerlst[3] = (z4,y4)
    stringerlst[4] = (z5,y5)
    stringerlst[5] = (z6,y6)
    stringerlst[6] = (z6,-y6)
    stringerlst[7] = (z5,-y5)
    stringerlst[8] = (z4,-y4)
    stringerlst[9] = (z3,-y3)
    stringerlst[10] = (z,-y)
    return stringerlst

def stringerylstgen():
    stringerylst = []
    stringerlst = stringerlstgen()
    for i in range(len(stringerlst)):
        stringerylst += [stringerlst[i][1]]
    return stringerylst

def stringerzlstgen():
    stringerzlst = []
    stringerlst = stringerlstgen()
    for i in range(len(stringerlst)):
        stringerzlst += [stringerlst[i][0]]
    return stringerzlst

def Astringer(w,h,t):
    #Input width, heights and thickness. Output the area of the stringers
    return((w+h)*t)

def Apanel(l,t):
    #input the length and thickness of the panels, output is the area of the panel
    return (l*t)

def Asemicircle(R,r):
    #input the radius and thickness of the semi circle, output the area of the semi circle
    return((np.pi/2)*(R**2-r**2))

def Aspar(h,t):
    #input the height and thickness of the spar. Output the area of the spar.
    return (h*t)

def Alstgen(Nstringer,Nspar,Nsemi,Npanel):
    #input the amount of each element, output a list of all the areas.
    Alst = []
    for i in range(Nstringer):
        Alst += [Astringer(wst,hst,tst)]
    for i in range(Nspar):
        Alst += [Aspar(h,tsp)]
    for i in range(Nsemi):
        Alst += [Asemicircle(R,r)]        
    for i in range(Npanel):
        Alst += [Apanel(l,tsk)]
    return Alst

def ycentroidlstgen():
    stringerylst = stringerylstgen()
    ycentroidlst = stringerylst
    ycentroidlst += [0.]
    ycentroidlst += [0.]
    ycentroidlst += [l*m.sin(theta)/2]
    ycentroidlst += [-l*m.sin(theta)/2]
    return ycentroidlst

print(l*m.cos(theta)/2)

print(l*m.cos(theta)/2)

def zcentroidlstgen():
    stringerzlst = stringerzlstgen()
    zcentroidlst = stringerzlst
    zcentroidlst += [0.]
    zcentroidlst += [2*R/np.pi]
    zcentroidlst += [-l*m.cos(theta)/2]
    zcentroidlst += [-l*m.cos(theta)/2]
    return(zcentroidlst)

#stringer, spar, semi, top, bottom
def zcentroidgen():
    zcentroidlst = zcentroidlstgen()
    Alst = Alstgen(Nstringer,Nspar,Nsemi,Npanel) 
    lst = []
    for i in range(len(Alst)):
        lst += [Alst[i]*zcentroidlst[i]]
    return(sum(lst)/sum(Alst))

def ycentroidgen():
    ycentroidlst = ycentroidlstgen()
    Alst = Alstgen(Nstringer,Nspar,Nsemi,Npanel)
    lst = []
    for i in range(len(Alst)):
        lst += [Alst[i]*ycentroidlst[i]]
        return(sum(lst)/sum(Alst))
       
def MoIrotrec(b,h,theta):
    return ((b * h)*(b**2 * m.cos(m.radians(90)-theta)**2 + h**2 * m.sin(m.radians(90)-theta)**2)/12)

def MoIsemi(R,r):
    #input the outer and inner radius of semi circle, output the moment of inertia about the center! NOT CENTROID!
    return((np.pi/8)*(R**4-r**4))
    
def MoIzgen():
    MoIz = 0
    MoIz += MoIrotrec(tsp,h,0) #Spar
    MoIz += MoIsemi(R,r) #Leading edge
    MoIz += 2*MoIrotrec(l,tsk,theta) # Panels
    return MoIz

def MoIygen():
    MoIy = 0
    MoIy += MoIrotrec(h,tsp,0) #Spar
    MoIy += MoIsemi(R,r) #Leading edge
    MoIy += 2*MoIrotrec(tsk,l,m.radians(90)-theta) # Panels
    return MoIy  
        
zcentroid = zcentroidgen()
ycentroid = ycentroidgen()

ycentroidlst = ycentroidlstgen()
zcentroidlst = zcentroidlstgen()

#print(zcentroidlst,ycentroidlst)
print('zbar = ',zcentroid,', ybar = ',ycentroid)
Alst = Alstgen(Nstringer,Nspar,Nsemi,Npanel)
#steinerterms

#print("Alst = ", Alst)

Steinerzlst = []
for i in range(len(ycentroidlst)):
    Steinerzlst += [Alst[i]*abs(ycentroid-ycentroidlst[i])**2]

Steinerylst = []
for i in range(len(zcentroidlst)):
    Steinerylst += [Alst[i]*abs(zcentroid-zcentroidlst[i])**2]
    
#print(Steinerzlst)
#print(Steinerylst)
    
print("Izz_st = ", sum(Steinerzlst))
print("Iyy_st = ", sum(Steinerylst))
print("Iyy = ", sum(Steinerylst)+MoIygen())
Steinerylst += [Asemicircle(R,r)*2*abs(zcentroid)*2*R/np.pi]
print(2*abs(zcentroid)*2*R/np.pi*Asemicircle(R,r))

#Izz = 1.028e-05
#Iyy = 5.377e-05
Izz = sum(Steinerzlst)+MoIzgen()
Iyy = sum(Steinerylst)+MoIygen()
print("Izz = ", sum(Steinerzlst)+MoIzgen())
print("Iyy = ", sum(Steinerylst)+MoIygen())

print(Iyy/5.377e-05)