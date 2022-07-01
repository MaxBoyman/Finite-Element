
import numpy as np

nodes_count=3
elements_count=2
A=7.65 #Surface
E=30e6	#elasticity coefficient
I=204   #second momentum of areas
L=np.array([10,9])*12 #Element lenght (UNIT : foot)
NNs=np.array([[1,2],[2,3]]) #Near Nodes
teta=np.array([0,270])    #Angle between element and global X
p=np.array([0,0,0,0,-4000,80000,0,0,0]) #force (X or Y)
bc=np.array([1,3])  #number of bc nodes


k_local=np.zeros((elements_count,6,6))
k_g=np.zeros((3*nodes_count,3*nodes_count))
k_g_bc=np.zeros((3*nodes_count,3*nodes_count))
u=np.zeros((3*nodes_count))
r=np.zeros((3*nodes_count))
def si(tet):
    s=np.sin(tet*np.pi/180)
    c=np.cos(tet*np.pi/180)
    T=np.array([[c,s,0,0,0,0],
                [-s,c,0,0,0,0],
                [0,0,1,0,0,0],
                [0,0,0,c,s,0],
                [0,0,0,-s,c,0],
                [0,0,0,0,0,1]])
    # print(T)
    return T

def base(A,E,L,I):
    B=np.array([[A*E/L,0,0,-(A*E/L),0,0],
                [0,12*((E*I)/L**3),6*((E*I)/L**2),0,-12*((E*I)/L**3),6*((E*I)/L**2)],
                [0,6*((E*I)/L**2),4*((E*I)/L),0,-6*((E*I)/L**2),2*((E*I)/L)],
                [-A*(E/L),0,0,A*(E/L),0,0],
                [0,-12*((E*I)/L**3),-6*((E*I)/L**2),0,12*((E*I)/L**3),-6*((E*I)/L**2)],
                [0,6*((E*I)/L**2),2*(E*I/L),0,-6*((E*I)/L**2),4*(E*I/L)]])
    # print(B)
    return B
for i in range(elements_count):
    k_local[i]=np.dot((np.dot(si(teta[i]).transpose(),base(A,E,L[i],I))),si(teta[i]))
# print(k_local,"\n")

for m in range(elements_count):
    i=NNs[m][0]
    j=NNs[m][1]
    k_g[3*i-3][3*i-3]=k_g[3*i-3][3*i-3]+k_local[m][0][0]
    k_g[3*i-3][3*i-2]=k_g[3*i-3][3*i-2]+k_local[m][0][1]
    k_g[3*i-3][3*i-1]=k_g[3*i-3][3*i-1]+k_local[m][0][2]
    k_g[3*i-3][3*j-3]=k_g[3*i-3][3*j-3]+k_local[m][0][3]
    k_g[3*i-3][3*j-2]=k_g[3*i-3][3*j-2]+k_local[m][0][4]
    k_g[3*i-3][3*j-1]=k_g[3*i-3][3*j-1]+k_local[m][0][5]

    k_g[3*i-2][3*i-3]=k_g[3*i-2][3*i-3]+k_local[m][1][0]
    k_g[3*i-2][3*i-2]=k_g[3*i-2][3*i-2]+k_local[m][1][1]
    k_g[3*i-2][3*i-1]=k_g[3*i-2][3*i-1]+k_local[m][1][2]
    k_g[3*i-2][3*j-3]=k_g[3*i-2][3*j-3]+k_local[m][1][3]
    k_g[3*i-2][3*j-2]=k_g[3*i-2][3*j-2]+k_local[m][1][4]
    k_g[3*i-2][3*j-1]=k_g[3*i-2][3*j-1]+k_local[m][1][5]

    k_g[3*i-1][3*i-3]=k_g[3*i-1][3*i-3]+k_local[m][2][0]
    k_g[3*i-1][3*i-2]=k_g[3*i-1][3*i-2]+k_local[m][2][1]
    k_g[3*i-1][3*i-1]=k_g[3*i-1][3*i-1]+k_local[m][2][2]
    k_g[3*i-1][3*j-3]=k_g[3*i-1][3*j-3]+k_local[m][2][3]
    k_g[3*i-1][3*j-2]=k_g[3*i-1][3*j-2]+k_local[m][2][4]
    k_g[3*i-1][3*j-1]=k_g[3*i-1][3*j-1]+k_local[m][2][5]

    k_g[3*j-3][3*i-3]=k_g[3*j-3][3*i-3]+k_local[m][3][0]
    k_g[3*j-3][3*i-2]=k_g[3*j-3][3*i-2]+k_local[m][3][1]
    k_g[3*j-3][3*i-1]=k_g[3*j-3][3*i-1]+k_local[m][3][2]
    k_g[3*j-3][3*j-3]=k_g[3*j-3][3*j-3]+k_local[m][3][3]
    k_g[3*j-3][3*j-2]=k_g[3*j-3][3*j-2]+k_local[m][3][4]
    k_g[3*j-3][3*j-1]=k_g[3*j-3][3*j-1]+k_local[m][3][5]

    k_g[3*j-2][3*i-3]=k_g[3*j-2][3*i-3]+k_local[m][4][0]
    k_g[3*j-2][3*i-2]=k_g[3*j-2][3*i-2]+k_local[m][4][1]
    k_g[3*j-2][3*i-1]=k_g[3*j-2][3*i-1]+k_local[m][4][2]
    k_g[3*j-2][3*j-3]=k_g[3*j-2][3*j-3]+k_local[m][4][3]
    k_g[3*j-2][3*j-2]=k_g[3*j-2][3*j-2]+k_local[m][4][4]
    k_g[3*j-2][3*j-1]=k_g[3*j-2][3*j-1]+k_local[m][4][5]

    k_g[3*j-1][3*i-3]=k_g[3*j-1][3*i-3]+k_local[m][5][0]
    k_g[3*j-1][3*i-2]=k_g[3*j-1][3*i-2]+k_local[m][5][1]
    k_g[3*j-1][3*i-1]=k_g[3*j-1][3*i-1]+k_local[m][5][2]
    k_g[3*j-1][3*j-3]=k_g[3*j-1][3*j-3]+k_local[m][5][3]
    k_g[3*j-1][3*j-2]=k_g[3*j-1][3*j-2]+k_local[m][5][4]
    k_g[3*j-1][3*j-1]=k_g[3*j-1][3*j-1]+k_local[m][5][5]
    
# print(k_g)

k_g_bc=k_g
for i in range(len(bc)):
    temp=np.array((2*nodes_count))
    k_g_bc[3*bc[i]-3]=0
    k_g_bc[3*bc[i]-2]=0
    k_g_bc[3*bc[i]-1]=0
    k_g_bc[3*bc[i]-3][3*bc[i]-3]=1
    k_g_bc[3*bc[i]-2][3*bc[i]-2]=1
    k_g_bc[3*bc[i]-1][3*bc[i]-1]=1


# print(k_g_bc)

u=np.linalg.solve(k_g_bc,p)

# print(u)

r=np.matmul(u,k_g)-p

print(r)