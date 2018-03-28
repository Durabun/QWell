import cmath
import numpy as np
import math
#from decimal import Decimal
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os, time, glob

#Defines constants in SI units
h = (6.626*10**-34)/(2*np.pi)

def dij(ki,kj):
        return 0.5*np.matrix(((1+kj/ki,1-kj/ki),(1-kj/ki,1+kj/ki)))
def Pi(ki,xi):
        return np.matrix(((np.exp(-(1j)*ki*xi),0),(0,np.exp((1j)*ki*xi))))
def Pj(kj,xj):
        return np.matrix(((np.exp((1j)*kj*xj),0),(0,np.exp(-(1j)*kj*xj))))

def K(m,E,V):
        return cmath.sqrt((2*m*(E-V)/h**2))

def CalculateCoeff(m,E,n,S,V,X):
	
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
	F = np.absolute(Mat[0,0]**-2)
	R = F*np.absolute(Mat[1,0]**2)
	return (F,R)
	#Prints the transmission and reflections coefficients
#	print(np.absolute(Mat[0,0]**-2))
#	F = np.absolute(Mat[0,0]**-2)
#	print(F*np.absolute(Mat[1,0]**2))

def LenBar(X):
	j = len(X)
	x = 0
	for i in range(j):
		x = x+ X[i]
		i = i +1
	return x

def plots(T,R,k,x):
	X = np.arange(-10*x*10**10,0,0.1)
	Z = np.arange(x*10**10,10*x*10**10,0.1)
	Y1 = np.exp((1j)*k*10**-10*X)
	Y2 = R*np.exp((-1j)*k*10**-10*X)
	Y3 = T*np.exp((1j)*k*10**-10*Z)
	plt.figure()
	axes = plt.gca()
	axes.set_xlim([-10*x*10**10,10*x*10**10])
	plt.plot(X,Y1,X,Y2,Z,Y3)
	plt.title("Multi-Barrier Wave Function")
	plt.xlabel('Angstrom')
	plt.ylabel('Probability Amplitude')
	plt.legend(["Incoming","Reflected","Transmitted"])
	if not os.path.isdir('static'):
       		os.mkdir('static')
	else:
        # Remove old plot files
        	for filename in glob.glob(os.path.join('static', '*.png')):
            		os.remove(filename)
    # Use time since Jan 1, 1970 in filename in order make
    # a unique filename that the browser has not chached
    	plotfile = os.path.join('static', str(time.time()) + '.png')
    	plt.savefig(plotfile)
    	return plotfile

