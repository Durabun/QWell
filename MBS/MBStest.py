import numpy as np
import scipy
import cmath
import math
from decimal import Decimal


#Test program for a single square barrier

#Defining parameters in NU

m = float((9.11*10**-31)/(1.8*10**-27)) #E
a = float((0.5*10**-9)/(0.197*10**(-15))) #1/E
E = float((0.75*1.602*10**-19)/(1.602*10**-10)) #E
V = float((1.602*10**-19)/(1.602*10**-10)) #E
k1 = cmath.sqrt((2*m*E)) #E
k2 = cmath.sqrt((2*m*(E-V))) #E

print(m)
print(k1)
print(k2)
print(k1/k2)
print(k2/k1)

def P12(k1,k2):
	return 0.5*np.matrix(((1+k2/k1,1-k2/k1),(1-k2/k1,1+k2/k1)))
def P21(k1,k2):
	return 0.5*np.matrix(((1+k1/k2,1-k1/k2),(1-k1/k2,1+k1/k2)))
def P2(k2,x):
	return np.matrix(((np.exp(-(1j)*k2*x),0),(0,np.exp((1j)*k2*x))))

def P1(k1,x):
	return np.matrix(((np.exp((1j)*k1*x),0),(0,np.exp(-(1j)*k1*x))))

A = P12(k1,k2)
B = P21(k1,k2)
C = P2(k2,a)
D = P1(k1,a)

print(A)

M = np.dot(np.dot(A,B,C),D)

print(M)

m1 = M[0,0]
M1 = cmath.polar(m1)
print(M1[0]) #Why does this only come out as 1???

X = cmath.polar((4*k1*k2)/(k1+k2)**2)

print(X)
