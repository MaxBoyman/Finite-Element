import numpy as np

nodes_count=4
elements_count=3
E=70e9 #Yang Module
Gama=5.4 #mass of lenght unit
L=np.array([0.1,0.1,0.1])   #Element Lenght
A=20e-4 #Surface


k_local=np.zeros((elements_count,nodes_count,nodes_count))
k_base=np.array([[1,-1],[-1,1]])
m_local=np.zeros((elements_count,nodes_count,nodes_count))
m_base=np.array([[2,1],[1,2]])
k_g=np.zeros((nodes_count,nodes_count))
k_g_bc=np.zeros((nodes_count,nodes_count))
m_g=np.zeros((nodes_count,nodes_count))
m_g_bc=np.zeros((nodes_count,nodes_count))


for i in range(elements_count):
        k_local[i]=np.lib.pad(((A*E/L[i])*k_base), ((i,elements_count-1-i),(i,elements_count-1-i)), 'constant', constant_values=(0))
        m_local[i]=np.lib.pad(((Gama*L[i]/6)*m_base), ((i,elements_count-1-i),(i,elements_count-1-i)), 'constant', constant_values=(0))


# print(k_local)
# print(f_local)

for i in range(0,elements_count,1):
    k_g=np.add(k_g,k_local[i])
    m_g=np.add(m_g,m_local[i])

# print(k_g)
# print(m_g)
k_g_bc=k_g
m_g_bc=m_g

k_g_bc=np.delete(k_g, 0,0)
m_g_bc=np.delete(m_g,0,0)
k_g_bc=np.delete(k_g_bc, 0,1)
m_g_bc=np.delete(m_g_bc,0,1)

# print(k_g_bc)
# print(m_g_bc)
m_g_bc_inv=np.linalg.inv(m_g_bc)

w=np.dot(m_g_bc_inv,k_g_bc)
# print(w)

W=np.linalg.eigvals(w)
print(np.sqrt(W))
