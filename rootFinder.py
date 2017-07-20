
import numpy as np
from os import system
from scipy.optimize import*


h = 6.626*10**(-34) #These are some constants
h_bar = h/(2*np.pi) #Specifically, Planck's Constant

#The fn below defines the even root function
def Evenroot(x,r):
	return x*((np.sin(x))/(np.cos(x)))-np.sqrt(r**2-x**2)


#The fn below defines the odd root function
def Oddroot(x,r):
	return -x*((np.cos(x))/(np.sin(x)))-np.sqrt(r**2-x**2)

#This fn is for determining the value of the even root curve. It is in fact non dimensionalized
def ErootY(x):
        return x*(np.sin(x))/(np.cos(x))

#This is the fn for determining the value of the odd root curve.
def OrootY(x):
        return -x*(np.cos(x))/(np.sin(x))

#This is to determine the propagation "K" vector. Units are in inverse Angstrom
def PropVec(x,a):
    return x/(a*10**(10))

#This is to determine the tunneling "k" vector (yeah, notation sucks). Units are in inverse Angstrom
def TunVec(y,a):
    return y/(a*10**(10))

#This function is for determining the energy of each state.
def Energy(y,p,V):
    return V*(y/p)**2

#The function below determines the radius of possible roots based on the parameters of the particle and the potential well
def radius(m,a,V):
	return np.sqrt((2*m*np.absolute(V)*a**2)/h_bar**2)

#Definition of the root radius. This will be determined by the user
R = radius((206.8*9.11*10**(-31)),(4.69*10**(-11)),(1.6022*10**(-19)))

i = 1
j = 0
while R > i*np.pi:
	i=i+1
while R > (1+2*j)*np.pi/2:
	j = j+1

print(str(i+j)+" root(s) in total.")


evenRootX = np.zeros(i)
evenRootY = np.zeros(i)
oddRootX = np.zeros(j)
oddRootY = np.zeros(j)
propagationVector = np.zeros(i+j)
tunnelVector = np.zeros(i+j)

ie = 0
jo = 0
ip = 0 #propagation vector index
it = 0 #tunnel vector index
dx = 0.1

#Test Cases for the initial case I worked on in my final
# In this case, there are 2 even roots and 1 odd root

def rootT(x):
	return x*((np.sin(x))/(np.cos(x)))-np.sqrt((4)**2-x**2)

def OrootT(x):
	return -x*((np.cos(x))/(np.sin(x)))-np.sqrt((4)**2-x**2)
evenRootX[ie] = newton(rootT,1+dx)
evenRootY[ie] = ErootY(evenRootX[ie])
propagationVector[ip] = PropVec(evenRootX[ie],(4.69*10**(-11)))
tunnelVector[it] = TunVec(evenRootY[ie],(4.69*10**(-11)))

#print newton(rootT,np.pi)
#print newton(OrootT,0.5*np.pi+3*dx)
#print newton(OrootT,1.5*np.pi+dx)

ie = 1
while i>ie:
	if R > ie*np.pi:
            evenRootX[ie] = newton(rootT,ie*np.pi+dx)
            evenRootY[ie] = ErootY(evenRootX[ie])
            ie = ie+1

while j>jo:
	if R > (1+2*jo)*np.pi/2:
            oddRootX[jo] = newton(OrootT,(1+2*jo)*np.pi/2+2*2*dx)
            oddRootY[jo] = OrootY(oddRootX[jo])
            jo=jo+1


print (evenRootX)
print (evenRootY)
print (oddRootX)
print (oddRootY)
print (propagationVector)
print(tunnelVector)
