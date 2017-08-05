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
while 1.1*np.pi > i*np.pi:
	i=i+1
while 1.1*np.pi > (1+2*j)*np.pi/2:
	j = j+1

print str(i+j)+" root(s) in total."


evenRoot = np.zeros(i)
oddRoot = np.zeros(j)

ie = 0
jo = 0
dx = 0.01
x0 = 0.01
#Test Cases for the initial case I worked on in my final
# In this case, there are 2 even roots and 1 odd root

#These are test root functions
#When I figure out how to make a function dynamic i.e. changing "r"
#Then I will go back to using the original root functions
def rootT(x):
	return x*((np.sin(x))/(np.cos(x)))-np.sqrt((1.1*np.pi)**2-x**2)

def OrootT(x):
	return -x*((np.cos(x))/(np.sin(x)))-np.sqrt((1.1*np.pi)**2-x**2)

while np.absolute(rootT(x0)) > 0.5:
	x0 = x0+dx
evenRoot[ie] = newton(rootT,x0)
#print newton(rootT,np.pi)
#print newton(OrootT,0.5*np.pi+3*dx)
#print newton(OrootT,1.5*np.pi+dx)

ie = 1
x0 = (1+2*jo)*np.pi/2
#while i>ie:
#	if R > ie*np.pi:
#		evenRoot[ie] = newton(rootT,ie*np.pi+dx)	
#		ie = ie+1
#while j>jo:
#	if R > (1+2*jo)*np.pi/2:
#		oddRoot[jo] = newton(OrootT,(1+2*jo)*np.pi/2+2*2*dx)
#		jo=jo+1

#This loop is for determing the X coords of the even and odd roots
#Which are then saved in appropriate arrays
while x0<1.1*np.pi:
	while np.absolute(OrootT(x0)) > 0.25:
		x0=x0+dx
	oddRoot[jo] = newton(OrootT,x0)
	jo=jo+1
	x0 = ie*np.pi
	print x0
	if x0<1.1*np.pi:
		while np.absolute(rootT(x0)) > 0.25:
			x0 = x0+dx
		evenRoot[ie] = newton(rootT,x0)
		ie = ie+1
		x0 = (1+2*jo)*np.pi/2
		print x0
	

print evenRoot
print oddRoot
print x0

