import numpy as np
from os import system
from scipy.optimize import*


h = 6.626*10**(-34) #These are some constants
h_bar = h/(2*np.pi) #Specifically, Planck's Constant

#The fn below defines the even root function
def Evenroot(r):
	return lambda x:x*((np.sin(x))/(np.cos(x)))-np.sqrt(r**2-x**2)


#The fn below defines the odd root function
def Oddroot(r):
	return lambda x:-x*((np.cos(x))/(np.sin(x)))-np.sqrt(r**2-x**2)

#Using a lambda function enables x as a variable while changing other variables

#The function below determines the radius of possible roots based on the parameters of the particle and the potential well
def radius(m,a,V):
	return np.sqrt((2*m*np.absolute(V)*a**2)/h_bar**2)
	
def energy(y,r,V):
	return np.absolute(V)*(y/r)**2

def radFunction(x,r):
	return np.sqrt(r**2-x**2)

def propVector(x,a):
	return (x/a)*10**(-10)

def tunnelVector(y,a):
	return (y/a)*10**(-10)

#Definition of the root radius. This will be determined by the user
R = radius((206.8*9.11*10**(-31)),(4.69*10**(-11)),(1.6022*10**(-19)))
x = 0
i = 1
j = 0
while R > i*np.pi:
	i=i+1
while R > (1+2*j)*np.pi/2:
	j = j+1

print (str(i+j)+" root(s) in total.")

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
#def rootT(x):
#	return x*((np.sin(x))/(np.cos(x)))-np.sqrt((R)**2-x**2)

#def OrootT(x):
#	return -x*((np.cos(x))/(np.sin(x)))-np.sqrt((R)**2-x**2)

	
ERoot = Evenroot(R)
ORoot = Oddroot(R)


while np.absolute(ERoot(x0)) > 0.5:
	x0 = x0+dx
evenRoot[ie] = newton(ERoot,x0)
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
while x0<R:
	while np.absolute(ORoot(x0)) > 0.25:
		x0=x0+dx
	oddRoot[jo] = newton(ORoot,x0)
	jo=jo+1
	x0 = ie*np.pi
#	print (x0)
	if x0<R:
		while np.absolute(ERoot(x0)) > 0.25:
			x0 = x0+dx
		evenRoot[ie] = newton(ERoot,x0)
		ie = ie+1
		x0 = (1+2*jo)*np.pi/2
#		print (x0)
	

print (evenRoot)
print (oddRoot)
#print (x0)

#initialize Even and Odd States
Even_State = np.zeros(i)
Odd_State = np.zeros(j)
#initialize Even and Odd tunnel vectors
tunnelEven = np.zeros(i)
tunnelOdd = np.zeros(j)
#initialize Even and Odd propagation vectors
propEven = np.zeros(i)
propOdd = np.zeros(j)


print (energy(radFunction(evenRoot[0],R),R,1))
print (tunnelVector(radFunction(evenRoot[0],R),4.69*10**(-11)))
print (propVector(evenRoot[0],4.69*10**(-11)))
#print (radFunction(evenRoot[0],R))