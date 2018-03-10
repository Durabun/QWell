import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import*
from scipy.special import*

#Defining constants by input
V = (1.602*10**(-19))*float(input("Depth? "))
m = (9.11*10**(-31))*float(input("Mass? "))
a = (10**(-10))*float(input("Half Width? "))
beta = V/a
hb = (6.626*10**(-34))/(2*np.pi)
alpha = ((2*m*beta)/(hb**2))**(1./3)

#Defining functions
#These can probably be tossed into a separate file

def SingAi(x):
	return airy(-x)[0]
def SingAip(x):
	return airy(-x)[1]
def SingBi(x):
	return airy(-x)[2]
def SingBip(x):
	return airy(-x)[3]

def DoubAi(alpha,a,x):
	return airy(alpha*a-x)[0]
def DoubAip(alpha,a,x):
	return airy(alpha*a-x)[1]
def DoubBi(alpha,a,x):
	return airy(alpha*a-x)[2]
def DoubBip(alpha,a,x):
	return airy(alpha*a-x)[3]

def Kvector(alpha,a,x):
	return (alpha**3*a-alpha**2*x)**(1./2)

#These are the root functions
def EvenR(alpha,a,x):
	return Kvector(alpha,a,x)*(SingAip(x)*DoubBi(alpha,a,x)-SingBip(x)*DoubAi(alpha,a,x))+alpha*(SingAip(x)*DoubBip(alpha,a,x)-SingBip(x)*DoubAip(alpha,a,x))

def OddR(alpha,a,x):
	return Kvector(alpha,a,x)*(SingAi(x)*DoubBi(alpha,a,x)-SingBi(x)*DoubAi(alpha,a,x))+alpha*(SingAi(x)*DoubBip(alpha,a,x)-SingBi(x)*DoubAip(alpha,a,x))

#Defining Energy and the rest of the constants
def Energy(alpha,b,x):
	return (x*b/alpha)

def TrueKV(alpha,V,m,E):
	return (2*m*(V-E)/hb**2)**(1./2)

def deltaE(alpha,beta,E):
	return -(airy(-alpha*E/beta)[1])/(airy(-alpha*E/beta)[3])

def deltaO(alpha,beta,E):
	return -(airy(-alpha*E/beta)[0])/(airy(-alpha*E/beta)[2])

def eta(alpha,beta,delta,a,E,K):
	return (airy(alpha*(a-E/beta))[0]+delta*airy(alpha*(a-E/beta))[2])*np.exp(K*a)

#Defining the Wave Functions
#Even
def Psi1E(eta,K):
	return lambda x: eta*np.exp(10**(-10)*K*x)

def Psi2E(alpha,beta,delta,E):
	return lambda x: airy(alpha*10**(-10)*(-x-10**(10)*E/beta))[0]+delta*airy(alpha*10**(-10)*(-x-10**(10)*E/beta))[2]

def Psi3E(alpha,beta,delta,E):
	return lambda x: airy(alpha*10**(-10)*(x-10**(10)*E/beta))[0]+delta*airy(alpha*10**(-10)*(x-10**(10)*E/beta))[2]

def Psi4E(eta,K):
	return lambda x: eta*np.exp(-10**(-10)*K*x)

#Odd
def Psi1O(eta,K):
	return lambda x: -eta*np.exp(10**(-10)*K*x)

def Psi2O(alpha,beta,delta,E):
	return lambda x: -airy(alpha*10**(-10)*(-x-10**(10)*E/beta))[0]-delta*airy(alpha*10**(-10)*(-x-10**(10)*E/beta))[2]

def Psi3O(alpha,beta,delta,E):
	return lambda x: airy(alpha*10**(-10)*(x-10**(10)*E/beta))[0]+delta*airy(alpha*10**(-10)*(x-10**(10)*E/beta))[2]

def Psi4O(eta,K):
	return lambda x: eta*np.exp(-10**(-10)*K*x)

#Root finding
S = 0
x0 = 0.5
dx = 0.01
T = 1

#This section iterates through the even root function
#When there is an error with the function it breaks through the loop
#This is because the bound root function just stops after a certain point
while (T==1):
	try:
		EvenR(alpha,a,x0+dx)
	except:
		T=0
		break
	#This uses a sign approach for determining roots
	#If the sign of EvenR(x0) is the same as the sign of EvenR(x0+dx)
	#Then x0 increases by dx. Eventually the signs will not be the same
	while (np.sign(EvenR(alpha,a,x0))==np.sign(EvenR(alpha,a,x0+dx))):
		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
		except:
			T=0
			break
		x0 = x0+dx
	T = 0
#There is almost always 1 even root. The above finds where it is located


#Same thing but for the odd roots
while (T==0):
	try:
		OddR(alpha,a,x0+dx)
		OddR(alpha,a,x0)
	except:
		T=1
		break
	while (np.sign(OddR(alpha,a,x0)) == np.sign(OddR(alpha,a,x0+dx))):
		x0 = x0+dx
		try:
			OddR(alpha,a,x0+dx)
			OddR(alpha,a,x0)
		except:
			T=1
			break
	#At the end of this loop, the number of roots S increases
	S = S +1
	try:
		EvenR(alpha,a,x0+dx)
		EvenR(alpha,a,x0)
	except:
		T=1
		break
	while(np.sign(EvenR(alpha,a,x0)) == np.sign(EvenR(alpha,a,x0+dx))):
		x0 = x0+dx
		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
		except:
			T=1
			break
	#This entire process determines how many roots there are (even and odd)
	S = S + 1

print("There are "+str(S)+" roots")
i = 0
x0 = 0.5

#These arrays hold the Energies and other constants for each state
#Based on the flow of the program, these can all be defined one fell swoop
Roots = np.zeros(S)
Energies = np.zeros(S)
Kvect = np.zeros(S)
Delta = np.zeros(S)
Eta = np.zeros(S)

#Thus begins another loop to go through the roots and calculate the energies at each root using Newtons method
while (i != S):
	if (i%2 == 0):
		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
		except:
			break
		while (np.sign(EvenR(alpha,a,x0)) == np.sign(EvenR(alpha,a,x0+dx))):
			try:
				EvenR(alpha,a,x0)
				EvenR(alpha,a,x0+dx)
			except:
				break
			x0 = x0+dx
			try:
				EvenR(alpha,a,x0)
				EvenR(alpha,a,x0+dx)
			except:
				break
		
		if (i==S):
			break
		#At this x0, the roots, energies, K vectors, and constants are calc'd
		Roots[i] = newton(lambda x: EvenR(alpha,a,x),x0)
		Energies[i] = Energy(alpha,beta,Roots[i])
		Kvect[i] = TrueKV(alpha,V,m,Energies[i])
		Delta[i] = deltaE(alpha,beta,Energies[i])
		Eta[i] = eta(alpha,beta,Delta[i],a,Energies[i],Kvect[i])

		i = i+1

	try:
		EvenR(alpha,a,x0)
		EvenR(alpha,a,x0+dx)
		OddR(alpha,a,x0)
		OddR(alpha,a,x0+dx)
	except:
		break

	if (i%2 != 0):
		try:
			OddR(alpha,a,x0)
			OddR(alpha,a,x0+dx)
		except:
			break
		while (np.sign(OddR(alpha,a,x0)) == np.sign(OddR(alpha,a,x0+dx))):
			try:
				OddR(alpha,a,x0)
				OddR(alpha,a,x0+dx)
			except:
				break
			x0 = x0+dx
			try:
				OddR(alpha,a,x0)
				OddR(alpha,a,x0+dx)
			except:
				break
		if (i == S):
			break
		#Same done here but for the odd roots
		Roots[i] = newton(lambda x: OddR(alpha,a,x),x0)
		Energies[i] = Energy(alpha,beta,Roots[i])
		Kvect[i] = TrueKV(alpha,V,m,Energies[i])
		Delta[i] = deltaO(alpha,beta,Energies[i])
		Eta[i] = eta(alpha,beta,Delta[i],a,Energies[i],Kvect[i])

		i = i+1

		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
			OddR(alpha,a,x0)
			OddR(alpha,a,x0+dx)
		except:
			break

#Displaying important numbers
print(Roots)
print(Energies)
print(Kvect)
print(Delta)
print(Eta)
	
#This section is identical to the finite well sim
#Z is used to iterate through the bound functions
#C is used to iterate through the colors for the plots
#dy is the function offset
z = 0
c = 0
dy = 0

r1 = np.arange((-2*a*10**(10)),(-a*10**(10)),0.00001)
r2 = np.arange((-a*10**(10)),0,0.00001)
r3 = np.arange(0,(a*10**(10)),0.00001)
r4 = np.arange((a*10**(10)),(2*a*10**(10)),0.00001)
color = ["r","b","g","c","m","y","k"]

#Plot function
while z < S:
	if c == 6:
		c = 0
	#Plotting odd functions
	if z%2 == 1:
		dy = (1/(1.602*10**(-19)))*Energies[z]
		plt1 = Psi1O(Eta[z],Kvect[z]) 
		plt2 = Psi2O(alpha,beta,Delta[z],Energies[z])
		plt3 = Psi3O(alpha,beta,Delta[z],Energies[z])
		plt4 = Psi4O(Eta[z],Kvect[z])

		plt.plot(r1,plt1(r1)+dy,color[c],r2,plt2(r2)+dy,color[c],r3,plt3(r3)+dy,color[c],r4,plt4(r4)+dy,color[c])
	#Plotting even functions
	else:
		dy = (1/(1.602*10**(-19)))*Energies[z]
		plt1 = Psi1E(Eta[z],Kvect[z]) 
		plt2 = Psi2E(alpha,beta,Delta[z],Energies[z])
		plt3 = Psi3E(alpha,beta,Delta[z],Energies[z])
		plt4 = Psi4E(Eta[z],Kvect[z])

		plt.plot(r1,plt1(r1)+dy,color[c],r2,plt2(r2)+dy,color[c],r3,plt3(r3)+dy,color[c],r4,plt4(r4)+dy,color[c])

	z = z+1
	c = c+1
#Presents all the plots
plt.show()
	












