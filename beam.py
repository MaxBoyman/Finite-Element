
import numpy as np

nodes_count=3
elements_count=2
A=6650e-6 #Surface
E=200e9	#elasticity coefficient
I=118.6e-6   #second momentum of areas
L=np.array([5,2.5])#*12 #Element lenght (UNIT : foot)
NNs=np.array([[1,2],[2,3]]) #Near Nodes
teta=np.array([0,0])    #Angle between element and global X
p=np.array([0,0,0,39062,-31250,13021]) #force (X or Y)
bc=np.array([1])  #number of bc nodes


k_local=np.zeros((elements_count,4,4))
k_g=np.zeros((2*nodes_count,2*nodes_count))
k_g_bc=np.zeros((2*nodes_count,2*nodes_count))
u=np.zeros((2*nodes_count))
r=np.zeros((2*nodes_count))
def si(tet):
    s=np.sin(tet*np.pi/180)
    c=np.cos(tet*np.pi/180)
    T=np.array([[c,s,0,0],
                [-s,c,0,0],
                [0,0,c,s],
                [0,0,-s,c]])
    # print(T)
    return T
def base(E,L,I):
    B=np.array([[12*((E*I)/L**3),6*((E*I)/L**2),-12*((E*I)/L**3),6*((E*I)/L**2)],
                [6*((E*I)/L**2),4*((E*I)/L),-6*((E*I)/L**2),2*((E*I)/L)],
                [-12*((E*I)/L**3),-6*((E*I)/L**2),12*((E*I)/L**3),-6*((E*I)/L**2)],
                [6*((E*I)/L**2),2*(E*I/L),-6*((E*I)/L**2),4*(E*I/L)]])
    # print(B)
    return B

for i in range(elements_count):
    k_local[i]=np.dot((np.dot(si(teta[i]).transpose(),base(E,L[i],I))),si(teta[i]))
# print(k_local,"\n")

for m in range(elements_count):
    i=NNs[m][0]
    j=NNs[m][1]
    k_g[2*i-2][2*i-2]=k_g[2*i-2][2*i-2]+k_local[m][0][0]
    k_g[2*i-2][2*i-1]=k_g[2*i-2][2*i-1]+k_local[m][0][1]
    k_g[2*i-2][2*j-2]=k_g[2*i-2][2*j-2]+k_local[m][0][2]
    k_g[2*i-2][2*j-1]=k_g[2*i-2][2*j-1]+k_local[m][0][3]

    k_g[2*i-1][2*i-2]=k_g[2*i-1][2*i-2]+k_local[m][1][0]
    k_g[2*i-1][2*i-1]=k_g[2*i-1][2*i-1]+k_local[m][1][1]
    k_g[2*i-1][2*j-2]=k_g[2*i-1][2*j-2]+k_local[m][1][2]
    k_g[2*i-1][2*j-1]=k_g[2*i-1][2*j-1]+k_local[m][1][3]

    k_g[2*j-2][2*i-2]=k_g[2*j-2][2*i-2]+k_local[m][2][0]
    k_g[2*j-2][2*i-1]=k_g[2*j-2][2*i-1]+k_local[m][2][1]
    k_g[2*j-2][2*j-2]=k_g[2*j-2][2*j-2]+k_local[m][2][2]
    k_g[2*j-2][2*j-1]=k_g[2*j-2][2*j-1]+k_local[m][2][3]

    k_g[2*j-1][2*i-2]=k_g[2*j-1][2*i-2]+k_local[m][3][0]
    k_g[2*j-1][2*i-1]=k_g[2*j-1][2*i-1]+k_local[m][3][1]
    k_g[2*j-1][2*j-2]=k_g[2*j-1][2*j-2]+k_local[m][3][2]
    k_g[2*j-1][2*j-1]=k_g[2*j-1][2*j-1]+k_local[m][3][3]
    
print(k_g)

k_g_bc=k_g
for i in range(len(bc)):
    temp=np.array((2*nodes_count))
    k_g_bc[2*bc[i]-2]=0
    k_g_bc[2*bc[i]-1]=0
    k_g_bc[2*bc[i]-2][2*bc[i]-2]=1
    k_g_bc[2*bc[i]-1][2*bc[i]-1]=1
    k_g_bc[2]=0
    k_g_bc[2][2]=1

print(k_g_bc)

u=np.linalg.solve(k_g_bc,p)

print(u.reshape(2*nodes_count,1))

r=np.matmul(u,k_g)-p

print(r)