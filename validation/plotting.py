import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d




########## INPUT

all_nodes = np.genfromtxt('input\\nodes_xyz.txt', delimiter =',', comments="*")
nodes_1_16 = np.genfromtxt('input\\nodes_1-16.txt', delimiter =',', comments="*")
nodes_skin = np.genfromtxt('input\\nodes_skin.txt', delimiter =',', comments="*")

all_elements = np.genfromtxt('input\\elements_4_nodes.txt', delimiter =',', comments="*")



x= np.zeros(len(all_nodes))
y= np.zeros(len(all_nodes))
z= np.zeros(len(all_nodes))

x2= np.zeros(len(nodes_1_16))
y2= np.zeros(len(nodes_1_16))
z2= np.zeros(len(nodes_1_16))

x3 = []
y3 = []
z3 = []

x_ele= np.zeros(len(all_elements))
y_ele= np.zeros(len(all_elements))
z_ele= np.zeros(len(all_elements))


i = 0
for row in all_nodes:
    x[i] = row[1]
    y[i] = row[2]
    z[i] = row[3]
    i +=1

j = 0
for row in nodes_1_16:
    x2[j] = row[1]
    y2[j] = row[2]
    z2[j] = row[3]
    j +=1

k = 0
for row in nodes_skin:
    for position in row:
        if position != 0:
            x3.append(x[int(position)-1])
            y3.append(y[int(position)-1])
            z3.append(z[int(position)-1])
            k +=1

l = 0
for row in all_elements:
    node_1 = int(row[1])
    node_2 = int(row[2])
    node_3 = int(row[3])
    node_4 = int(row[4])
    x_ele[l] = (all_nodes[node_1-1][1]+all_nodes[node_2-1][1]+all_nodes[node_3-1][1]+all_nodes[node_4-1][1])/4
    y_ele[l] = (all_nodes[node_1-1][2]+all_nodes[node_2-1][2]+all_nodes[node_3-1][2]+all_nodes[node_4-1][2])/4
    z_ele[l] = (all_nodes[node_1-1][3]+all_nodes[node_2-1][3]+all_nodes[node_3-1][3]+all_nodes[node_4-1][3])/4
    l +=1










##########  OUTPUT

#Bending-case van mises stress and S12 (Shear stress)
stress_b = np.genfromtxt('output\\b_ele_mises_s12.txt', skip_header= 3)
stress_b = stress_b[stress_b[:, 0].argsort()]

stress_b_2 = np.zeros((len(stress_b),3))
i=0
for row in stress_b:
    stress_b_2[i][0]=int(row[0])
    stress_b_2[i][1]=(row[2]+row[3])/2
    stress_b_2[i][2]=(row[4]+row[5])/2
    i=i+1

#Jammed bending case van mises stress and S12 (Shear stress)
stress_jb = np.genfromtxt('output\\jb_ele_mises_s12.txt', skip_header= 3)
stress_jb = stress_jb[stress_jb[:, 0].argsort()]

stress_jb_2 = np.zeros((len(stress_jb),3))
i=0
for row in stress_jb:
    stress_jb_2[i][0]=int(row[0])
    stress_jb_2[i][1]=(row[2]+row[3])/2
    stress_jb_2[i][2]=(row[4]+row[5])/2
    i=i+1
    
#Jammed straight case van mises stress and S12 (Shear stress)
stress_js = np.genfromtxt('output\\js_ele_mises_s12.txt', skip_header= 3)
stress_js = stress_js[stress_js[:, 0].argsort()]

stress_js_2 = np.zeros((len(stress_js),3))
i=0
for row in stress_js:
    stress_js_2[i][0]=int(row[0])
    stress_js_2[i][1]=(row[2]+row[3])/2
    stress_js_2[i][2]=(row[4]+row[5])/2
    i=i+1



#Bending-case van mises stress and S12 (Shear stress)
deflec_b = np.genfromtxt('output\\b_node_deflec.txt', skip_header= 3)




# Finds the element with the highest von mises stress and highest Shear stress S12
max_stress = np.amax(stress_jb_2, axis=0)
max_stress_ele = np.flip(np.array(np.where(stress_jb_2==np.amax(stress_jb_2, axis=0))))

# Mises Stress in cross-section


cut_x_mises = x_ele[max_stress_ele[1][1]]


print("Max. von Mises stress in element: ", max_stress_ele[1][1]+1)
print("with von Mises stress of: ", max_stress[1])
print("at x= ", x_ele[max_stress_ele[1][1]])
print("at y= ", y_ele[max_stress_ele[1][1]])
print("at z= ", z_ele[max_stress_ele[1][1]])
print()



y_c_sec_mises = np.array([])
z_c_sec_mises = np.array([])
mises_jb_c_sec = np.array([])




i=0
for x in x_ele:
    if x==cut_x_mises:
        y_c_sec_mises = np.append(y_c_sec_mises, y_ele[i])
        z_c_sec_mises = np.append(z_c_sec_mises, z_ele[i])
        mises_jb_c_sec = np.append(mises_jb_c_sec, stress_jb_2[i, 1])
    i=i+1

min_mises_c_sec = np.amin(mises_jb_c_sec)
index_min_mises_c_sec = int(np.where(mises_jb_c_sec==min_mises_c_sec)[0])

print("Min. von Mises stress in the same cross section of: ", min_mises_c_sec)
print("at x= ", cut_x_mises)
print("at y= ", y_c_sec_mises[index_min_mises_c_sec])
print("at z= ", z_c_sec_mises[index_min_mises_c_sec])
print()
print()
print()




# Shear in cross-section

cut_x_s12 = x_ele[max_stress_ele[1][2]]

print("Max. Shear stress in element: ", max_stress_ele[1][2]+1)
print("with Shear stress of: ", max_stress[2])
print("at x= ", x_ele[max_stress_ele[1][2]])
print("at y= ", y_ele[max_stress_ele[1][2]])
print("at z= ", z_ele[max_stress_ele[1][2]])
print()


y_c_sec_s12 = np.array([])
z_c_sec_s12 = np.array([])
s12_jb_c_sec = np.array([])


i=0
for x in x_ele:
    if x==cut_x_s12:
        y_c_sec_s12 = np.append(y_c_sec_s12, y_ele[i])
        z_c_sec_s12 = np.append(z_c_sec_s12, z_ele[i])
        s12_jb_c_sec = np.append(s12_jb_c_sec, stress_jb_2[i, 2])
    i=i+1
        

min_s12_c_sec = np.amin(s12_jb_c_sec)
index_min_s12_c_sec = int(np.where(s12_jb_c_sec==min_s12_c_sec)[0])

print("Min. Shear stress in the same cross section of: ", min_s12_c_sec)
print("at x= ", cut_x_s12)
print("at y= ", y_c_sec_s12[index_min_s12_c_sec])
print("at z= ", z_c_sec_s12[index_min_s12_c_sec])






######## PLOTTING

fig1 = plt.figure()
ax1 = plt.axes()
ax1.invert_xaxis()
img = plt.scatter(z_c_sec_mises, y_c_sec_mises, c = mises_jb_c_sec, cmap ='RdYlGn_r')

fig1.colorbar(img)
ax1.set_title("Von Mises stresses, Case: jammed bending")
ax1.set_xlabel('z')
ax1.set_ylabel('y')



fig2 = plt.figure()
ax2 = plt.axes()
ax2.invert_xaxis()
img = plt.scatter(z_c_sec_mises, y_c_sec_mises, c = s12_jb_c_sec, cmap ='RdYlGn_r')

fig2.colorbar(img)
ax2.set_title("S12 Shear stress, Case: jammed bending")
ax2.set_xlabel('z')
ax2.set_ylabel('y')




fig3 = plt.figure()
ax3 = plt.axes()
ax3.invert_xaxis()
img = plt.scatter(z_c_sec_mises, y_c_sec_mises, c = s12_jb_c_sec, cmap ='RdYlGn_r')

fig3.colorbar(img)
ax3.set_title("Deflection, Case: jammed bending")
ax3.set_xlabel('z')
ax3.set_ylabel('y')


#cmap = plt.get_cmap("jet")
#img = plt.plot(z_mises_c_sec, y_mises_c_sec)
#m = cm.ScalarMappable(cmap=plt.get_cmap("jet"))
#img = plt.scatter(z_mises_c_sec, y_mises_c_sec, c = stress_c_sec)
#fig.colorbar(img)

#img = ax.scatter3D(x, y, z)
#img2 = ax.scatter(x2, y2, z2)
#img3 = ax.scatter(x3, y3, z3)
#img4 = ax.scatter3D(x_ele, y_ele, z_ele)

#ax.set_zlabel('z')
plt.show()