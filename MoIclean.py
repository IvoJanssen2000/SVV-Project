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

#Computing handy parameters
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

def Astringer(w,h,t):
    #Input width, heights and thickness. Output the area of the stringers
    return((w+h)*t)

def Arect(b,h):
    #input the width and length of a rectangle, output the area of the rectangle
    return (b*h )

def Asemicircle(R,r):
    #input the outer and inner radius of the semi circle, output the area of the semi circle
    return((np.pi/2)*(R**2-r**2))

def MoIrotrec(b,h,theta):
    #input the width, height and rotation angle of a rectangle, output the moment of inertia about its horizontal axis.
    return ((b * h)*(b**2 * m.cos(m.radians(90)-theta)**2 + h**2 * m.sin(m.radians(90)-theta)**2)/12)

def MoIsemi(R,r):
    #input the outer and inner radius of a semi circle, output the moment of inertia about the center.
    return((np.pi/8)*(R**4-r**4))

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

def MoIzgen(): #creates Izz
    MoIz = 0
    MoIz += MoIrotrec(tsp,h,0) #Spar
    MoIz += MoIsemi(R,r) #Leading edge
    MoIz += 2*MoIrotrec(l,tsk,theta) # Panels
    return MoIz

def MoIygen(): #creates Iyy
    MoIy = 0
    MoIy += MoIrotrec(h,tsp,0) #Spar
    MoIy += MoIsemi(R,r) #Leading edge
    MoIy += 2*MoIrotrec(tsk,l,m.radians(90)-theta) # Panels
    return MoIy  

Steinerzlst = [] #creates a list with the steiner terms for Izz
for i in range(len(ycentroidlst)):
    Steinerzlst += [Alst[i]*abs(ycentroid-ycentroidlst[i])**2]
Steinerylst = [] #creates a list with the steiner terms for Iyy
for i in range(len(zcentroidlst)):
    Steinerylst += [Alst[i]*abs(zcentroid-zcentroidlst[i])**2]

Izz = sum(Steinerzlst)+MoIzgen() #sum the Izz steiner terms and the Izz contributions of all the elements
Iyy = sum(Steinerylst)+MoIygen() #sum the Iyy steiner terms and the Iyy     contributions of all the elements
print("Izz = ", Izz)
print("Iyy = ", Iyy)