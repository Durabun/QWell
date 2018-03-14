import cmath
import numpy as np
import math
from decimal import Decimal


#Test program for multiple square barriers
#No inputs (yet)
#So the user will need to change the following constants


m = float((9.11*10**-31)) 
a = float(10**(-10))
b = a
c = a
E = float(0.9613*10**-19)
V1 = float(1.2817*10**-19) 
#V2 = float(1.602*10**-19)
V2 = E
V3 = 1*V1
h = (6.626*10**-34)/(2*np.pi)
#k1 = cmath.sqrt((2*m*E/h**2)) 
#k2 = cmath.sqrt((2*m*(E-V)/h**2))



#Defining the transfer matrices and K vector
def dij(ki,kj):
	return 0.5*np.matrix(((1+kj/ki,1-kj/ki),(1-kj/ki,1+kj/ki)))
def Pi(ki,xi):
	return np.matrix(((np.exp(-(1j)*ki*xi),0),(0,np.exp((1j)*ki*xi))))
def Pj(kj,xj):
	return np.matrix(((np.exp((1j)*kj*xj),0),(0,np.exp(-(1j)*kj*xj))))

def K(m,E,V):
	return cmath.sqrt((2*m*(E-V)/h**2))

#Outside Identity matrix
Q = np.matrix(((1,0),(0,1)))
#Number of wells
n = 3
#Number of Matrices
s = 2*(n+1)
#Array of potential energies and thicknesses
#These need to be changed manually right now
V = np.array([0,V1,V2,V3])
X = np.array([a,b,c])

#Defining indeces. For example, matrix d12
#Would correspond to the continuity between waves 1 and 2
#i would be 1 and j would be 2 in this case
i=0
j=1

#Defines the transfer matrix which is dependent on the amount of barriers
for i in range(s):
#	print(Q)
	if i == n:
		j = 0
		N = dij(K(m,E,V[i]),K(m,E,V[j]))
		#These are the last 2 matrices which "return" to the original wave
		N = np.dot(N,Pj(K(m,E,V[j]),X[j+1]))
		break
	#Progresses through the first s-2 matrices
	M = dij(K(m,E,V[i]),K(m,E,V[j]))
	M = np.dot(M,Pi(K(m,E,V[j]),X[i]))
	Q = np.dot(Q,M)
	i = j
	j = j+1

#Transfer matrix
Mat = np.dot(Q,N)


print(Mat)

#Prints the transmission and reflections coefficients
print(np.absolute(Mat[0,0]**-2))
F = np.absolute(Mat[0,0]**-2)
print(F*np.absolute(Mat[1,0]**2))

















