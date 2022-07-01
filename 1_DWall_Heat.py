import numpy as np

nodes_count=4
elements_count=3
T_Base=200  #Base Temperature
T_Fluid=30  #Fluid Temp
K=np.array([0.08,0.074,0.72])   #Thermal Conductivity
h=40    #Convection Coefficient
L=np.array([5e-2,15e-2,10e-2])
A=1   #surface area
p=0   #mohit
NNs=np.array([[1,2],[2,3],[3,4]]) #Near Nodes
bc=np.array([1,4])  #number of bc nodes


k_local=np.zeros((elements_count,nodes_count,nodes_count))
k1_base=np.array([[1,-1],[-1,1]])
k2_base=np.array([[2,1],[1,2]])
k_g=np.zeros((nodes_count,nodes_count))
k_g_bc=np.zeros((nodes_count,nodes_count))
t=np.zeros((nodes_count))
f_base=np.array([1,1]).reshape(2,1)
f_local=np.zeros((elements_count,nodes_count,1))
f_g=np.zeros((nodes_count,1))

for i in range(elements_count):
    if i<elements_count-1:
        k_local[i]=np.lib.pad(((K[i]*A/L[i])*k1_base)+((h*p*L[i]/6)*k2_base), ((i,elements_count-1-i),(i,elements_count-1-i)), 'constant', constant_values=(0))
        f_local[i]=np.lib.pad(((h*p*L[i]*T_Fluid/2)*f_base),  ((i,2-i),(0,0)), 'constant', constant_values=(0))
    else:
        k_local[i]=np.lib.pad(((K[i]*A/L[i])*k1_base)+((h*p*L[i]/6)*k2_base)+np.array([[0,0],[0,h*A]]), ((i,elements_count-1-i),(i,elements_count-1-i)), 'constant', constant_values=(0))
        f_local[i]=np.lib.pad(((h*p*L[i]*T_Fluid/2)*f_base)+np.array([0,(h*A*T_Fluid)]).reshape(2,1),  ((i,elements_count-1-i),(0,0)), 'constant', constant_values=(0))

# print(k_local)
# print(f_local)

for i in range(0,elements_count,1):
    k_g=np.add(k_g,k_local[i])
    f_g=np.add(f_g,f_local[i])

print(k_g)
print(f_g)
k_g[0][0]=1
k_g[0][1]=0
f_g[0]=T_Base

t=np.linalg.solve(k_g,f_g)

print(t)