
import numpy as np

nodes_count=4
elements_count=4
A=0.002 #Surface
E=70e9	#elasticity coefficient
L=np.array([1.5,np.sqrt((1.5**2)+1),1.5,1])#*12 #Element lenght (UNIT : foot)
NNs=np.array([[1,2],[2,3],[3,4],[4,2]]) #Near Nodes
teta=np.array([0,(np.arctan(1.5)+np.pi/2)*180/np.pi,0,90])    #Angle between element and global X
p=np.array([0,0,0,0,0,0,0,-10000]) #force (X or Y)
bc=np.array([1,3])  #number of bc nodes


k_local=np.zeros((elements_count,4,4))
k_g=np.zeros((2*nodes_count,2*nodes_count))
k_g_bc=np.zeros((2*nodes_count,2*nodes_count))
u=np.zeros((2*nodes_count))
r=np.zeros((2*nodes_count))
def si(tet):
    s=np.sin(tet*np.pi/180)
    c=np.cos(tet*np.pi/180)
    T=np.array([[c**2,c*s,-(c**2),-(s*c)],
                [s*c,s**2,-(s*c),-(s**2)],
                [-(c**2),-(c*s),c**2,c*s],
                [-(s*c),-(s**2),s*c,s**2]])
    # print(T)
    return T

for i in range(elements_count):
    k_local[i]=((A*E)/L[i])*si(teta[i])
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


print(k_g_bc)

u=np.linalg.solve(k_g_bc,p)

print(u.reshape(2*nodes_count,1))

r=np.matmul(u,k_g)-p

print(r)