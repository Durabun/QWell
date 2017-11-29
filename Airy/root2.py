import numpy as np
import matplotlib.pyplot as plt
from scipy.special import*
from scipy.optimize import*

r = np.arange(-10,0,0.1)
X = np.arange(0.1,10,0.1)

V = 1.602*10**(-19)
a = 0.5*10**(-10)
beta = V/a
m = 206*9.11*10**(-31)
hb = (6.626*10**(-34))/(2*np.pi)
alpha = ((2*m*beta)/hb**2)**(1./3)

def Sing():
	return lambda x: airy(-x)

def Doub(alpha,a):
	return lambda x: airy(alpha*a-x)

def KVec(alpha,a):
	return lambda x: (alpha**3*a-alpha**2*x)**(1./2)

A = Sing()
B = Doub(alpha,a)
K = KVec(alpha,a)

def SingAi(x):
	return airy(-x)[0]
def SingBi(x):
	return airy(-x)[2]
def SingAip(x):
	return airy(-x)[1]
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

def EvenR(alpha,a,x):
	return Kvector(alpha,a,x)*(SingAip(x)*DoubBi(alpha,a,x)-SingBip(x)*DoubAi(alpha,a,x))+alpha*(SingAip(x)*DoubBip(alpha,a,x)-SingBip(x)*DoubAip(alpha,a,x))

def OddR(alpha,a,x):
	return Kvector(alpha,a,x)*(SingAi(x)*DoubBi(alpha,a,x)-SingBi(x)*DoubAi(alpha,a,x))+alpha*(SingAi(x)*DoubBip(alpha,a,x)-SingBi(x)*DoubAip(alpha,a,x))

def Energy(alpha,b,x):
	return x*b/alpha

def TrueKV(alpha,V,m,E):
	return (2*m*(V-E)/hb**2)**(1./2)


#These are the Even Constants
def deltaE(alpha,beta,E):
	return -(airy(-alpha*E/beta)[1])/(airy(-alpha*E/beta)[3])

def etaE(alpha,beta,delta,a,E,K):
	return (airy(alpha*(a-E/beta))[0]+delta*airy(alpha*(a-E/beta))[2])*np.exp(K*a)


#These are the Odd Constants
def deltaO(alpha,beta,E):
	return -(airy(-alpha*E/beta)[0])/(airy(-alpha*E/beta)[2])
#Eta is the same for Even and Odd functions
def etaO(alpha,beta,delta,a,E,K):
	return (airy(alpha*(a-E/beta))[0]+float(delta)*airy(alpha*(a-E/beta))[2])*np.exp(K*a)


#Even Wavefunction
def Psi1E(eta,K):
	return lambda x: eta*np.exp(K*x)

def Psi2E(alpha,beta,delta,E):
	return lambda x: airy(alpha*(-x-E/beta))[0]+delta*airy(alpha*(-x-E/beta))[2]

def Psi3E(alpha,beta,delta,E):
	return lambda x: airy(alpha*(x-E/beta))[0]+delta*airy(alpha*(x-E/beta))[2]

def Psi4E(eta,K):
	return lambda x: eta*np.exp(-K*x)


#Odd Wavefunction
def Psi1O(eta,K):
	return lambda x: -eta*np.exp(K*x)

def Psi2O(alpha,beta,delta,E):
	return lambda x: -airy(alpha*(-x-E/beta))[0]-delta*airy(alpha*(-x-E/beta))[2]

def Psi3O(alpha,beta,delta,E):
	return lambda x: airy(alpha*(x-E/beta))[0]+delta*airy(alpha*(x-E/beta))[2]

def Psi4O(eta,K):
	return lambda x: eta*np.exp(-K*x)


S = 0
x0 = 0.5
dx = 0.01
T = 1

while (T==1):
	try:
		EvenR(alpha,a,x0+dx)
	except:
		T=0
		break
	while (np.sign(EvenR(alpha,a,x0))==np.sign(EvenR(alpha,a,x0+dx))):
		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
		except:
			T=0
			break
		x0 = x0+dx
	T = 0

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
	S = S + 1

print("There are "+str(S)+" roots")
i = 0
x0 = 0.5
Roots = np.zeros(S)
Energies = np.zeros(S)
Kvect = np.zeros(S)
Delta = np.zeros(S)
Eta = np.zeros(S)

while (i != S):
	if (i%2 == 0):
		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
		except:
			print("Fail1")
			break
		while (np.sign(EvenR(alpha,a,x0))==np.sign(EvenR(alpha,a,x0+dx))):
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

		Roots[i] = newton(lambda x: EvenR(alpha,a,x),x0)
		Energies[i] = Energy(alpha,beta,Roots[i]) 
		Kvect[i] = TrueKV(alpha,V,m,Energies[i])
		Delta[i] = deltaE(alpha,beta,Energies[i])
		Eta[i] = etaE(alpha,beta,Delta[i],a,Energies[i],Kvect[i])
		
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
			OddR(alpha,a,x0+dx)
			OddR(alpha,a,x0)
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
		if(i == S):
			break

		Roots[i] = newton(lambda x: OddR(alpha,a,x),x0)
		Energies[i] = Energy(alpha,beta,Roots[i])
		Kvect[i] = TrueKV(alpha,V,m,Energies[i])
		Delta[i] = deltaO(alpha,beta,Energies[i])
		Eta[i] = etaO(alpha,beta,Delta[i],a,Energies[i],Kvect[i])

		i = i+1

		try:
			EvenR(alpha,a,x0)
			EvenR(alpha,a,x0+dx)
			OddR(alpha,a,x0)
			OddR(alpha,a,x0+dx)
		except:
			break

print(Roots)
print(Energies)
print(Kvect)
print(Delta)
print(Eta)

Even = K(X)*(A(X)[1]*B(X)[2]-A(X)[3]*B(X)[0])+alpha*(A(X)[1]*B(X)[3]-A(X)[3]*B(X)[1])

Odd = K(X)*(A(X)[0]*B(X)[2]-A(X)[2]*B(X)[0])+alpha*(A(X)[0]*B(X)[3]-A(X)[2]*B(X)[1])

#r1 = np.arange(-2*a,-a,0.0001)
#r2 = np.arange(-a,0)
#r3 = np.arange(0,a,0.0001)
#r4 = np.arange(a,2*a)

#Atry = 0.5

#Test = np.arange(a*10**(10),2*a*10**(10),0.0001)

#plt.plot(X,Even)
#plt.plot(X,Odd)
#plt.plot((0,5),(0,0))
#E1 = Psi1E(Eta[0],Kvect[0])
#E2 = Psi2E(alpha,beta,Delta[0],Energies[0])
#E3 = Psi3E(alpha,beta,Delta[0],Energies[0])
#E4 = Psi4E(Eta[0],Kvect[0])



#plt.plot(Test,Eta[0]*np.exp())
#plt.show()
