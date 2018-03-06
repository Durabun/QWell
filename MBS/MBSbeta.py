import cmath
import numpy as np
import math
from decimal import Decimal


#Test program for a single square barrier


m = float((9.11*10**-31)) 
a = float(10**(-10))
E = float(0.8*1.602*10**-19) 
V = float(1.602*10**-19) 
h = (6.626*10**-34)/(2*np.pi)
#k1 = cmath.sqrt((2*m*E/h**2)) 
#k2 = cmath.sqrt((2*m*(E-V)/h**2))




def dij(ki,kj):
	return 0.5*np.matrix(((1+kj/ki,1-kj/ki),(1-kj/ki,1+kj/ki)))
def Pi(ki,xi):
	return np.matrix(((np.exp(-(1j)*ki*xi),0),(0,np.exp((1j)*ki*xi))))
def Pj(kj,xj):
	return np.matrix(((np.exp((1j)*kj*xj),0),(0,np.exp(-(1j)*kj*xj))))

def K(m,E,V):
	return cmath.sqrt((2*m*(E-V)/h**2))

n = 2
s = 2*(n+1)
V = np.array((1.2*10**-19,1.6*10**-19))
thickness = np.array((10**-10,10**-10))

kVect = np.zeros(n+1,dtype=np.complex)
k=cmath.sqrt((2*m*E/h**2))
kVect[0] = k

i=0
j=1

for i in range(n):
	kVect[j]=K(m,E,V[j-1])
	j=j+1
		
i=0
j=1

for i in range(n):
	if i == n:
		j = 0
		D = dij(k[i],k[j])
		Y = D
		P = Pj(k[j],thickness[i])
		Y = np.dot(D,P)
		break
	D = dij(k[i],k[j])
	X = d
	P = Pj(k[j],thickness[i])
	X = np.dot(X,P)
	i = i+1
	j = j+1


























