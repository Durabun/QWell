import cmath
import numpy as np
import math
from decimal import Decimal

#Defines constants in SI units
h = (6.626*10**-34)/(2*np.pi)
m = float(input("mass? ")*9.11*10**-31)
E = float(input("energy? ")*1.602*10**-19)

#Defines the number of barriers and boundaries
n = input("Number of barriers? ")
S = 2*(n+1)

#initializes heights and thicknesses as fillable arrays
V = np.zeros(n+1)
X = np.zeros(n)

#Loop iterates through the arrays to fill in the heights and thicknesses
l = 1
for l in range(n):
	#The first spot in V corresponds to the free particle, so V0 = 0
	V[l+1] = float(input("Height of barrier "+str(l+1)+": ")*1.602*10**-19)
	X[l] = float(input("Thickness of barrier "+str(l+1)+": ")*10**-10)
	l=l+1

#Defining the transfer matrices and K vector
def dij(ki,kj):
        return 0.5*np.matrix(((1+kj/ki,1-kj/ki),(1-kj/ki,1+kj/ki)))
def Pi(ki,xi):
        return np.matrix(((np.exp(-(1j)*ki*xi),0),(0,np.exp((1j)*ki*xi))))
def Pj(kj,xj):
        return np.matrix(((np.exp((1j)*kj*xj),0),(0,np.exp(-(1j)*kj*xj))))

def K(m,E,V):
        return cmath.sqrt((2*m*(E-V)/h**2))

#Defining Identity matrix outside
Q = np.matrix(((1,0),(0,1)))

i = 0
j = 1

#Defines the transfer matrix which is dependent on the amount of barriers
for i in range(S):
#       print(Q)
        if i == n:
                j = 0
                N = dij(K(m,E,V[i]),K(m,E,V[j]))
                #These are the last 2 matrices which "return" to the original wave
                N = np.dot(N,Pj(K(m,E,V[j]),X[j]))
                break
        #Progresses through the first s-2 matrices
        M = dij(K(m,E,V[i]),K(m,E,V[j]))
        M = np.dot(M,Pi(K(m,E,V[j]),X[i]))
        Q = np.dot(Q,M)
        i = j
        j = j+1

#Transfer Matrix
Mat = np.dot(Q,N)

#Prints the transmission and reflections coefficients
print(np.absolute(Mat[0,0]**-2))
F = np.absolute(Mat[0,0]**-2)
print(F*np.absolute(Mat[1,0]**2))


