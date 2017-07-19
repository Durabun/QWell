
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

print str(i+j)+" root(s) in total."


evenRoot = np.zeros(i)
oddRoot = np.zeros(j)

ie = 0
jo = 0
dx = 0.1

#Test Cases for the initial case I worked on in my final
# In this case, there are 2 even roots and 1 odd root

def rootT(x):
	return x*((np.sin(x))/(np.cos(x)))-np.sqrt((4)**2-x**2)

def OrootT(x):
	return -x*((np.cos(x))/(np.sin(x)))-np.sqrt((4)**2-x**2)
evenRoot[ie] = newton(rootT,1+dx)
#print newton(rootT,np.pi)
#print newton(OrootT,0.5*np.pi+3*dx)
#print newton(OrootT,1.5*np.pi+dx)

ie = 1
while i>ie:
	if R > ie*np.pi:
		evenRoot[ie] = newton(rootT,ie*np.pi+dx)	
		ie = ie+1
while j>jo:
	if R > (1+2*jo)*np.pi/2:
		oddRoot[jo] = newton(OrootT,(1+2*jo)*np.pi/2+2*2*dx)
		jo=jo+1


print evenRoot
print oddRoot

