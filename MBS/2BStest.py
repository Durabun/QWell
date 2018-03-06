import numpy as np
import scipy
import cmath
import math
from decimal import Decimal

#Test for a double barrier system
#Note
#a,b,c,... should probably be the thicknesses of the barriers instead of the locations since I redefined some of the coordinates of the barriers

m = float((9.11*10**-31)) 
a = float(10**(-10))
b = float(10**-10)
E = float(0.9613*10**-19) 
V1 = float(1.2817*10**-19) 
V2 = float(1.602*10**-19)
h = (6.626*10**-34)/(2*np.pi)
k1 = cmath.sqrt((2*m*E/h**2)) 
k2 = cmath.sqrt((2*m*(E-V1)/h**2)) 
k3 = cmath.sqrt((2*m*(E-V2)/h**2)) 

def dij(ki,kj):
	return 0.5*np.matrix(((1+kj/ki,1-kj/ki),(1-kj/ki,1+kj/ki)))
def Pi(ki,xi):
	return np.matrix(((np.exp(-(1j)*ki*xi),0),(0,np.exp((1j)*ki*xi))))
def Pj(kj,xj):
	return np.matrix(((np.exp((1j)*kj*xj),0),(0,np.exp(-(1j)*kj*xj))))

d12 = dij(k1,k2)
P2 = Pi(k2,a)
d23 = dij(k2,k3)
P3 = Pi(k3,b)
d31 = dij(k3,k1)
P1 = Pj(k1,b)


print(d12)
print(P2)
print(d23)
print(P3)
print(d31)
print(P1)


A = np.dot(d12,P2)
B = np.dot(d23,P3)
C = np.dot(d31,P1)

D = np.dot(A,B)
M = np.dot(D,C)

print(M)

print(np.absolute(M[0,0]**-2))
F = np.absolute(M[0,0]**-2)
print(F*np.absolute(M[1,0]**2))

