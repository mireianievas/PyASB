import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Se trata de ajustar X = A*T
# valores de x -
#t = np.random.normal(3,15,20)
#u = np.random.normal(2,12,20)

t = [np.random.random()*10 for s in range(50)]
u = [np.random.random()*10 for s in range(50)]


T = np.array([t,u,np.sin(u)])

# A = [3 0 4; 1 1 0]
x = 3*T[0]+4*T[2]
y = T[0] + T[2]

# Generamos ruido
max_noise = 1e-4
ruido1 = np.random.normal(0,max_noise,len(t))
ruido2 = np.random.normal(0,max_noise,len(u))

x = x + ruido1
y = y + ruido2

X = np.array([x,y]).transpose()

n = len(T)
k = len(T[1])

# modelo

M = np.ones(n)
for m in range(1,k):
	#print M
	#print T[:,m]
	#M = np.append([M],[T[:,m]],axis=0)
	M = np.column_stack([M,T[:,m]])

M = M.transpose()

#print M
#print np.dot(M.transpose(),M)
#print M.transpose()

b = np.linalg.inv(np.dot(M.transpose(),M))
b = np.dot(b,M.transpose())
b = np.dot(b,X)

#b = np.dot(M.transpose(),X)
#b = np.dot(np.linalg.inv(np.dot(M.transpose(),M)),b)

SSE = (np.dot(M,b)-X)
SSE = np.dot(SSE.transpose(),SSE)

#print SSE

Xmedia = np.mean(X)
SSR = (np.dot(M,b)-Xmedia)
SSR = np.dot(SSR.transpose(),SSR)

#print SSR


SST = SSE + SSR

v1 = k
v2 = n-k-1

s2 = SSE/(n-k-1)
F = (SSR/k)/s2

print 'Curva ajustada: '
print b.transpose()











