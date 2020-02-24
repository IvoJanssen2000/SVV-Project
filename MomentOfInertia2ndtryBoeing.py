# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:46:19 2020

@author: marti
"""

import math as m
import numpy as np

Nstringer = 0
Nspar = 1
Npanel = 2
Nsemi = 1
Ca = 0.605      #m
h = 0.205       #m 
tsk = 0.0011    #m
tsp = 0.0028    #m
tst = 0.0012    #m
hst = 0.016     #m
wst = 0.019      #m

#Computing handy parameters
R = h/2
theta = m.atan(R/(Ca-R))
l = np.sqrt((Ca-R)**2 + R**2)
L = np.pi*R + 2* l
tnew = tsk+15*((wst+hst)*tst)/L
print(tnew)
r = h/2 - tnew

def Astringer(w,h,t):
    #Input width, heights and thickness. Output the area of the stringers
    return((w+h)*t)

def Arect(b,h):
    #input the width and length of a rectangle, output the area of the rectangle
    return (b*h )

def Asemicircle(R,r):
    #input the outer and inner radius of the semi circle, output the area of the semi circle
    return((np.pi/2)*(R**2-r**2))

def Alstgen(Nstringer,Nspar,Nsemi,Npanel):
    #input the amount of each element, output a list of all the areas.
    Alst = []
    for i in range(Nstringer):
        Alst += [Astringer(wst,hst,tst)]
    for i in range(Nspar):
        Alst += [Arect(h,tsp)]
    for i in range(Nsemi):
        Alst += [Asemicircle(R,r)]        
    for i in range(Npanel):
        Alst += [Arect(l,tsk)]
    return Alst

def MoIrotrec(b,h,theta):
    #input the width, height and rotation angle of a rectangle, output the moment of inertia about its horizontal axis.
    return ((b * h)*(b**2 * m.cos(m.radians(90)-theta)**2 + h**2 * m.sin(m.radians(90)-theta)**2)/12)

def MoIsemi(R,r):
    #input the outer and inner radius of a semi circle, output the moment of inertia about the center.
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

zcentroidlst = [ 0. , 0.06525353 , -0.25125 , -0.25125 ]
ycentroidlst = [ 0. , 0. , 0.05125 , -0.05125 ]

zcentroid  = -0.13810173573534712
ycentroid = 0
Alst = Alstgen(Nstringer,Nspar,Nsemi,Npanel)
Steinerzlst = []
for i in range(len(ycentroidlst)):
    Steinerzlst += [Alst[i]*abs(ycentroid-ycentroidlst[i])**2]
Steinerylst = []
for i in range(len(zcentroidlst)):
    Steinerylst += [Alst[i]*abs(zcentroid-zcentroidlst[i])**2]

print("Izz = ", sum(Steinerzlst)+MoIzgen())
print("Iyy = ", sum(Steinerylst)+MoIygen())

Izz = sum(Steinerzlst)+MoIzgen()
Iyy = sum(Steinerylst)+MoIygen()
'Izz = 8.552855743273587e-06'
'Iyy = 4.968481864787321e-05'