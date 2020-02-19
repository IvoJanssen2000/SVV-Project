# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:35:11 2020

@author: marti
"""
import numpy as np
import math as m

#"Inputs:
#    Stringer locations
#    Stringer areas
#    Spar location
#    Spar area
#    Top surface area
#    Top surface location
#    Bottom surface area
#    Bottom surface location
#    Semi circle location
#    Semi circle area
#    
#    All scalars
#
#Output:
#    Moment of Inertia

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

#stringer centroid is 2.5 mm from base.
stcent = 0.0025
R = h/2
r = R - tsk
theta = m.atan(R/(Ca-R))
theta1 = m.atan(R/(Ca-R))
l = np.sqrt((Ca-R)**2 + R**2)
L = np.pi*R + 2* l
d = L/11
circ = 2*np.pi*R
angle = 360 * d/circ

#print(stiffangle)
theta2 = m.radians(angle)
theta3 = m.radians(90)
z = R*m.cos(m.radians(angle))
y = R*m.sin(m.radians(angle))
#print(z)
#print(y)

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

def Alstgen(Nstringer,Nspar,Nseml,Npanel):
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



Alst = Alstgen(Nstringer,Nspar,Nsemi,Npanel)
#print(Alst)

def MoIrect(b,h):
    #input the width and length of a rectangle, output the moment of inertia about its centroid.
    return((1/12)*b*h**3)
    
def MoIsemi(R,r):
    #input the outer and inner radius of semi circle, output the moment of inertia about the center! NOT CENTROID!
    return((np.pi/8)*(R**4-r**4))

def MoIstringerx(tst,hst,wst):
    #input the dimensions of the stringer and receive the complete MoIxx about its centroid.
    MoI1 = MoIrect(tst,hst)
    MoI2 = MoIrect(wst,tst)
    d1 = hst/2 - stcent
    d2 = stcent
    St1 = tst*hst*d1**2
    St2 = wst*tst*d2**2
    return(MoI1+MoI2+St1+St2)

def MoIstringery(tst,hst,wst):
    #input the dimensions of the stringer, output the MoIyy about its centroid.
    return(MoIrect(hst,tst) + MoIrect(tst,wst))

def MoIrotrec(b,h,theta):
    return ((b * h)*(b**2 * m.cos(m.radians(90)-theta)**2 + h**2 * m.sin(m.radians(90)-theta)**2)/12)

def MoIrotrec2(b,h,theta):
    return ((b**3 * h * m.sin(theta)**2 ) / 12)

print(MoIrotrec(l,tsk,-theta1))
print(MoIrotrec2(l,tsk,-theta1))


def MoIconv(Ixx,Iyy,theta):
    "This program calculates the Moment of Inertia when rotating the reference frame."
    "Input the moment of inertia about axis system x-y, and the rotation angle theta."
    "The output is the moment of inertia about the rotated axis system u-v."
    Iuu = (Ixx+Iyy)/2 + m.cos(2*theta)*(Ixx-Iyy)/2
    Ivv = (Ixx+Iyy)/2 - m.cos(2*theta)*(Ixx-Iyy)/2
    return Iuu,Ivv

def MoIlstgen(Nstringer,Nspar,Nsemi,Npanel,theta1,theta2,theta3):
    #input the amount of elements and an angle theta, output a list of the MoIxx and MoIyy of the different elements.
    MoIxlst = []
    MoIylst = []
    for i in range(Nstringer-3):
        MoIxlst += [0]
        MoIylst += [0]
        #MoIxlst += [MoIconv(MoIstringerx(tst,hst,wst),MoIstringery(tst,hst,wst),theta1)[0]]
        #MoIylst += [MoIconv(MoIstringerx(tst,hst,wst),MoIstringery(tst,hst,wst),theta1)[1]]
    for i in range(2):
        MoIxlst += [0]
        MoIylst += [0]
        #MoIxlst += [MoIconv(MoIstringerx(tst,hst,wst),MoIstringery(tst,hst,wst),theta2)[0]]
        #MoIylst += [MoIconv(MoIstringerx(tst,hst,wst),MoIstringery(tst,hst,wst),theta2)[1]]
    for i in range(1):
        MoIxlst += [0]
        MoIylst += [0]
        #MoIxlst += [MoIconv(MoIstringerx(tst,hst,wst),MoIstringery(tst,hst,wst),theta3)[0]]
        #MoIylst += [MoIconv(MoIstringerx(tst,hst,wst),MoIstringery(tst,hst,wst),theta3)[1]]        
    for i in range(Nspar):
        MoIxlst += [MoIrect(tsp,h)] 
        MoIylst += [MoIrect(h,tsp)]
    for i in range(Nsemi):
        MoIxlst += [MoIsemi(R,r)]  
        MoIylst += [MoIsemi(R,r)]
    for i in range(Npanel):
        MoIxlst += [MoIrotrec(l,tsk,theta)]
        MoIylst += [MoIrotrec(tsk,l,m.radians(90)-theta)]
        #MoIxlst += [MoIconv(MoIrect(l,tsk),MoIrect(tsk,l),theta1)[0]]
        #MoIylst += [MoIconv(MoIrect(l,tsk),MoIrect(tsk,l),theta1)[1]]

    return MoIxlst, MoIylst

MoIxlst,MoIylst = MoIlstgen(Nstringer,Nspar,Nsemi, Npanel,theta,theta2,theta3)
print(MoIxlst)
print(MoIconv(MoIrect(l,tsk),MoIrect(tsk,l),theta1)[0])


(zc,yc) = (-0.08327784945088734,0)

#print(z)
#print(y)
# order is: stringers, spar, semi, top, bottom.

zcentroidlst = [ 0.1215,      0.07675071 ,-0.0232709 , -0.12811921 ,-0.23296753, -0.33781584,
 -0.33781584, -0.23296753 ,-0.12811921 ,-0.0232709  , 0.07675071 , 0.,
  0., -0.1955     ,-0.1955    ]

ycentroidlst = [ 0.         , 0.09418648  ,0.11399895 , 0.08074782 , 0.04749669  ,0.01424556,
 -0.01424556 ,-0.04749669 ,-0.08074782 ,-0.11399895, -0.09418648  ,0.,
  0.,0.062 ,     -0.062     ]

def Steinerlstgen():
    Alst = Alstgen(Nstringer,Nspar,Nsemi,Npanel)
    Steinerxlst = []
    Steinerylst = []
    for i in range(Nstringer+Nspar+Nsemi+Npanel):
        #print(abs(z-zcentroidlst[i]))
        Steinerxlst += [Alst[i]*(ycentroidlst[i])**2]
        Steinerylst += [Alst[i]*(abs(zc-zcentroidlst[i]))**2]
    return Steinerxlst, Steinerylst

Steinerxlst,Steinerylst = Steinerlstgen()

#print(Steinerxlst,Steinerylst)

MoIzz = sum(MoIxlst) + sum(Steinerxlst)
MoIyy = sum(MoIylst) + sum(Steinerylst)


print("Izz = ",MoIzz)
print("Iyy = ",MoIyy)

#Izz = 1.028e-05
#Iyy = 8.651e-05

print(sum(Steinerxlst))
print(sum(Steinerylst))

testMoIxx = sum(Steinerxlst) + MoIrect(tsp,h) + MoIsemi(R,r) + 2*MoIrect(l,tsk)
print(testMoIxx)

