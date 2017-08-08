import numpy as np
from os import system
from scipy.optimize import*
import matplotlib.pyplot as plt
import quantumFn

h = 6.626*10**(-34) #These are some constants
h_bar = h/(2*np.pi) #Specifically, Planck's Constant

#This initializes the parameters of the problem
m =(9.11*10**(-31))*float(input("Mass of particle (X times the mass of electron): ")) #The mass is based on the electron
V = (1.602*10**(-19))*float(input("Depth of well (eV): "))
a = (10**(-10))*float(input("Half width of well centered at x=0 (A): "))

#Defines the "radius" based on the parameters
#Which will be used to determine the roots
R = quantumFn.radius(m,a,V)
x = 0
i = 1
j = 0

while R > i*np.pi:
	i = i+1
while R > (1+2*j)*np.pi/2:
	j = j+1
S = i+j

print (str(S)+" root(s) in total.")

#Initializes the "x" coord of the even and odd roots
evenRoot = np.zeros(i)
oddRoot = np.zeros(j)
#Used for indexing for the roots
ie = 0
jo = 0
#This is an interative root finder
#With an x0 used as an initial guess for the Newton function
dx = 0.01
x0 = 0.01

#Used the Even/Odd root functions to establish a function
#With the radius as the parameter and "x" as a variable
ERoot = quantumFn.Evenroot(R)
ORoot = quantumFn.Oddroot(R)

#Since there will always be one root, this while loop
#Searches for the first even root, and then puts it into
#The evenRoot array
while np.absolute(ERoot(x0)) > 0.5:
	x0 = x0 + dx
evenRoot[ie] = newton(ERoot,x0)

#If there is another even root, it will go into index 1
#For the even root array
#The next possible root will be when the radius is at least
#pi/2, so that is why x0 jumps to that value
ie = 1
x0 = (1+2*jo)*np.pi/2

#This loop will search for the remaining possible roots
#It goes in the order of looking for an odd root
#And then an even root.
while x0 < R:
	while np.absolute(ORoot(x0)) > 0.25:
		x0 = x0 + dx
	oddRoot[jo] = newton(ORoot,x0)
	jo = jo +1
	x0 = ie*np.pi #Resets the initial guess to next possible value
	if x0 < R:
		while np.absolute(ERoot(x0)) > 0.25:
			x0 = x0 +dx
		evenRoot[ie] = newton(ERoot,x0)
		ie = ie +1
		x0 = (1+2*jo)*np.pi/2 #Resets again

print (evenRoot)
print (oddRoot)

#Initialize the propagation and tunneling vectors
#At this point I consolidate the even and odd roots
#Into single arrays instead of having arrays for 
#Even and odd roots
propVector = np.zeros(S)
tunnVector = np.zeros(S)

i = 0
j = 1

#This is to define the vectors of each function
#THIS STILL DOES NOT WORK CORRECTLY
#IT MESSES UP WHEN THE NUMBER OF SOLUTIONS INCREASES
while i < evenRoot.size:
	propVector[2*i] = (10**(-10))*evenRoot[i]/(a)
	tunnVector[2*i] = (10**(-10))*(quantumFn.radFunction(evenRoot[i],R))/a
	i = i +1
while j-1 < oddRoot.size:
	propVector[j] = (10**(-10))*oddRoot[j-1]/(a)
	tunnVector[j] = ((10**(-10))*quantumFn.radFunction(oddRoot[j-1],R))/a
	j = j + 2

print ("Propagation 'k' vectors")
print (propVector)
print ("Tunneling 'K' vectors")
print (tunnVector)

#What needs to be done:
#1. Correct the vector definition
#2. Incorporate the definition of the B coefficients
#3. Incorporate the rest of index_plot.py


#Bcoeff = np.zeros(S)

Bcoeff = np.array([0.6318,0.617,0.4823])
#l = 0
#def B_Function(k,K)
#while l < S:
#	Bcoeff[l] = B_Function(prop[l],tunn[l])
#	l = l +1

print (Bcoeff)

z = 0
r1 = np.arange(-5,(-a*10**(10)),0.001)
r2 = np.arange((-a*10**(10)),(a*10**(10)),0.001)
r3 = np.arange((a*10**(10)),5,0.001)
color = ["r","b","g","c","m","y","k"] #This is the color array

while z < S:
	if z%2 == 1:
		plt1 = quantumFn.of1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
		plt2 = quantumFn.of2(Bcoeff[z],propVector[z])
		plt3 = quantumFn.of3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))

		plt.plot(r1,plt1(r1),color[z],r2,plt2(r2),color[z],r3,plt3(r3),color[z])

	else:
		plt1 = quantumFn.ef1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
		plt2 = quantumFn.ef2(Bcoeff[z],propVector[z])
		plt3 = quantumFn.ef3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))

		plt.plot(r1,plt1(r1),color[z],r2,plt2(r2),color[z],r3,plt3(r3),color[z])
	z = z + 1

plt.show()	
		
		







