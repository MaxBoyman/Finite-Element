import numpy as np

nodes_count=7
elements_count=6
u=np.array([0.775,0.525,0.34,0.24,0.16,0.11])  #Viscosity
ro=np.array([1254,1252,1249,1245,1241,1236])  #Density
K=np.array([0.08,0.074,0.72])   #Thermal Conductivity
w=40e-2    #channel width
L=np.array([1.5e-2,1.5e-2,1.5e-2,1.5e-2,1.5e-2,1.5e-2])
pd=-120  #pressure Drop dp/dx
NNs=np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7]]) #Near Nodes
bc=np.array([1,7])  #number of bc nodes


k_local=np.zeros((elements_count,nodes_count,nodes_count))
k_base=np.array([[1,-1],[-1,1]])
k_g=np.zeros((nodes_count,nodes_count))
k_g_bc=np.zeros((nodes_count,nodes_count))
v=np.zeros((nodes_count))   #speed
mdot=np.zeros((elements_count))   #fluids flow
f_base=np.array([1,1]).reshape(2,1)
f_local=np.zeros((elements_count,nodes_count,1))
f_g=np.zeros((nodes_count,1))
f_g_bc=np.zeros((nodes_count,1))

for i in range(elements_count):
        k_local[i]=np.lib.pad(((u[i]/L[i])*k_base), ((i,elements_count-1-i),(i,elements_count-1-i)), 'constant', constant_values=(0))
        f_local[i]=np.lib.pad(((-(pd*L[i]/2))*f_base),  ((i,elements_count-1-i),(0,0)), 'constant', constant_values=(0))


# print(k_local)
# print(f_local)

for i in range(0,elements_count,1):
    k_g=np.add(k_g,k_local[i])
    f_g=np.add(f_g,f_local[i])

# print(k_g)
# print(f_g)
k_g_bc=k_g
f_g_bc=f_g

k_g_bc[0][0]=1
k_g_bc[0][1]=0
k_g_bc[elements_count][elements_count]=1
k_g_bc[elements_count][elements_count-1]=0
f_g_bc[0]=0
f_g_bc[elements_count]=0

# print(k_g_bc)
# print(f_g_bc)

v=np.linalg.solve(k_g,f_g)

print("u:",v)

for i in range(elements_count):
    mdot[i]=(ro[i]*w*L[i])*((v[i]+v[i+1])/2)

print("M.",mdot)