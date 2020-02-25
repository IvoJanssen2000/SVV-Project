# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 10:53:36 2020

@author: elraw
"""
# len of region one = 5778
# len of region two = 855

import matplotlib.pyplot as plt
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D

# inut file
nodes    = np.genfromtxt('B737_input.INP', skip_header= 9 , skip_footer = 7996, delimiter = ','  )  
elements = np.genfromtxt('B737_input.INP', skip_header= 6598 , skip_footer = 1361, delimiter = ',' , dtype = int )  
     
# output file 
stress1_1    = np.genfromtxt('B737.RPT', skip_header= 20 , skip_footer = 53992) #len of region 1 : 0 -> 53992. len = 5778 +

stress1_2   = np.genfromtxt('B737.RPT', skip_header= 5778+20+18 , skip_footer = 53992-867)
stress= np.vstack((stress1_1,stress1_2))
stress_sorted = stress[stress[:,0].argsort()]  

#x= np.zeros((len(nodes)))
#y= np.zeros((len(nodes)))
#z= np.zeros((len(nodes)))


#for i in range(len(nodes)):
#    x[i] = nodes[i][1]
#    y[i] = nodes[i][2]
#    z[i] = nodes[i][3]

f = np.arange(0,len(nodes),1)

x= np.zeros((len(elements)))
y= np.zeros((len(elements)))
z= np.zeros((len(elements)))

i = 0
for row in elements:
    this_node = int(row [1]-1)
    x[i] = nodes[this_node][1]
    y[i] = nodes[this_node][2]
    z[i] = nodes[this_node][3]
    i +=1
    

segma_1 = stress_sorted[:,2]
segma_2 = stress_sorted[:,3]
#segma = 
fig = plt.figure()
ax = plt.axes(projection='3d')
img = ax.scatter(x,y,z, c = segma )
fig.colorbar(img)
plt.show()


#for i in range (len(nodes)):
 #   ele = np.where(elements == i)
 #   ele1 = ele[0]








''' 

stress_sorted = np.zeros((len(stress),6))
for i in range (len(stress)):
    index = np.where(stress[:,0] == i)
    print(index)
    lst.append(index) 
    #row = index[0][0]
    #stress_sorted[i,:] = stress[row]
    '''
    
    