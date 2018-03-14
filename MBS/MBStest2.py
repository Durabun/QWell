import numpy as np
import scipy
import cmath
import math
from decimal import Decimal


#Test program for a single square barrier


m = float((9.11*10**-31)) 
a = float(10**(-10))
E = float(0.81*1.602*10**-19) 
V = float(0.8*1.602*10**-19) 
h = (6.626*10**-34)/(2*np.pi)
k1 = cmath.sqrt((2*m*E/h**2)) 
k2 = cmath.sqrt((2*m*(E-V)/h**2))

print(k1)
print(k2)

def d12(k1,k2):
	return 0.5*np.matrix(((1+k2/k1,1-k2/k1),(1-k2/k1,1+k2/k1)))
def d21(k1,k2):
	#return 0.5*np.matrix(((1+(k1/k2),1-(k1/k2)),(1-(k1/k2),1+(k1/k2))))
	#return 0.5*np.matrix(((1+k2/k1,1-k2/k1),(1-k2/k1,1+k2/k1)))
	#return 0.5*(1+k1/k2)
	return 0.5*np.matrix(((1+k1/k2,1-k1/k2),(1-k1/k2,1+k1/k2)))

def P2(k2,x):
	return np.matrix(((np.exp(-(1j)*k2*x),0),(0,np.exp((1j)*k2*x))))
def P1(k1,x):
	return np.matrix(((np.exp((1j)*k1*x),0),(0,np.exp(-(1j)*k1*x))))


A = d12(k1,k2)
C = d21(k1,k2)
B = P2(k2,a)
D = P1(k1,a)

H = np.dot(A,B)
J = np.dot(C,D)

X = np.dot(H,J)

eheh = np.array((H,J))

FF = np.dot(eheh[0],eheh[1])

print(FF[0,0])
print(np.absolute(FF[0,0]**-2))
F = np.absolute(FF[0,0]**-2)

print(F*np.absolute(FF[1,0])**2)
