#QWell 0.5, the quantum well simulator!
#Made by Jonathan Eugenio
#And Mark Romero

#Please refer to the license in the github repo (https://github.com/Durabun/QWell.git)

import numpy as np
from os import system
from scipy.optimize import*
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import quantumFn
from matplotlib import animation

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

dx = 0.0001 #This value has to be smaller to resolve bigger
x0 = 0.000  #Scenarios!

#Used the Even/Odd root functions to establish a function
#With the radius as the parameter and "x" as a variable
ERoot = quantumFn.Evenroot(R)
ORoot = quantumFn.Oddroot(R)

#Since there will always be one root, this while loop
#Searches for the first even root, and then puts it into
#The evenRoot array
while np.absolute(ERoot(x0)) > 0.05: #This value has to
	x0 = x0 + dx		     #Be small to resolve
print (x0)			     #Smaller scenarios!
print (R)
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

#Initialize the propagation and tunneling vectors
#At this point I consolidate the even and odd roots
#Into single arrays instead of having arrays for 
#Even and odd roots
propVector = np.zeros(S)
tunnVector = np.zeros(S)

ev = 0
od = 1
Roots = np.zeros(S)
#This is to define the vectors of each function
while ev < S:
	Roots[ev] = evenRoot[int(ev/2)]
	ev = ev +2

while od < S:
	Roots[od] = oddRoot[int(od/2)]
	od = od + 2

print (Roots)

ev = 0
od = 1

while ev < S:
	propVector[ev] = (10**(-10))*Roots[ev]/a
	tunnVector[ev] = (10**(-10))*(quantumFn.radFunction(Roots[ev],R))/a
	ev = ev + 2

while od < S:
	propVector[od] = (10**(-10))*Roots[od]/a
	tunnVector[od] = (10**(-10))*(quantumFn.radFunction(Roots[od],R))/a
	od = od + 2

print ("Propagation 'k' vectors")
print (propVector)
print ("Tunneling 'K' vectors")
print (tunnVector)

Bcoeff = np.zeros(S)
StateEnergy = np.zeros(S)

l = 0
while l < S:
	Bcoeff[l] = quantumFn.Beven(a*10**(10),propVector[l],tunnVector[l])
	StateEnergy[l] = -1*quantumFn.energy(quantumFn.radFunction(Roots[l],R),R,V)
	l = l + 2
l = 1
while l < S:
	Bcoeff[l] = quantumFn.Bodd(a*10**(10),propVector[l],tunnVector[l])
	StateEnergy[l] = -1*quantumFn.energy(quantumFn.radFunction(Roots[l],R),R,V)
	l = l + 2

print ("Coefficients")
print (Bcoeff)
print ("Energies")
print (StateEnergy)

z = 0
r1 = np.arange((-5*a*10**(10)),(-a*10**(10)),0.001)
r2 = np.arange((-a*10**(10)),(a*10**(10)),0.001)
r3 = np.arange((a*10**(10)),(5*a*10**(10)),0.001)
color = ["r","b","g","c","m","y","k"] #This is the color array

c = 0
while z < S:
	if c == 6:
		c = 0
	if z%2 == 1:
		plt1 = quantumFn.of1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
		plt2 = quantumFn.of2(Bcoeff[z],propVector[z])
		plt3 = quantumFn.of3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
		
		plt.plot(r1,plt1(r1)+StateEnergy[z],color[c],r2,plt2(r2)+StateEnergy[z],color[c],r3,plt3(r3)+StateEnergy[z],color[c])
		
		

	else:
		plt1 = quantumFn.ef1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
		plt2 = quantumFn.ef2(Bcoeff[z],propVector[z])
		plt3 = quantumFn.ef3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))

		plt.plot(r1,plt1(r1)+StateEnergy[z],color[c],r2,plt2(r2)+StateEnergy[z],color[c],r3,plt3(r3)+StateEnergy[z],color[c])
		
		
	z = z + 1
	c = c + 1

plt.xlabel("Distance (in Angstrom)")
plt.show()
