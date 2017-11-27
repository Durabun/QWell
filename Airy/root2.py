import numpy as np
import matplotlib.pyplot as plt
from scipy.special import*
from scipy.optimize import*

r = np.arange(-10,0,0.1)
X = np.arange(0.1,10,0.1)
x = 0
S = 0

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

def Kvector(x):
	return (alpha**3*a-alpha**2*x)**(1./2)

def EvenR(alpha,a,x):
	return Kvector(x)*(SingAip(x)*DoubBi(alpha,a,x)-SingBip(x)*DoubAi(alpha,a,x))+alpha*(SingAip(x)*DoubBip(alpha,a,x)-SingBip(x)*DoubAip(alpha,a,x))

def OddR(alpha,a,x):
	return Kvector(x)*(SingAi(x)*DoubBi(alpha,a,x)-SingBi(x)*DoubAi(alpha,a,x))+alpha*(SingAi(x)*DoubBip(alpha,a,x)-SingBi(x)*DoubAip(alpha,a,x))

#print("This is a test Eval")
#print(np.sign(EvenR(alpha,a,1)))
#print(np.sign(EvenR(alpha,a,1.1)))

xE = 0.5
xO = 0.5
dx = 0.1

while (np.sign(EvenR(alpha,a,xE))==np.sign(EvenR(alpha,a,xE+dx))):
	xE = xE+dx

print(xE)

while (np.sign(OddR(alpha,a,xO))==np.sign(OddR(alpha,a,xO+dx))):
	xO = xO+dx

print(xO)



print(newton(lambda x:EvenR(alpha,a,x),xE))

print(newton(lambda x:OddR(alpha,a,x),xO))

Even = K(X)*(A(X)[1]*B(X)[2]-A(X)[3]*B(X)[0])+alpha*(A(X)[1]*B(X)[3]-A(X)[3]*B(X)[1])

Odd = K(X)*(A(X)[0]*B(X)[2]-A(X)[2]*B(X)[0])+alpha*(A(X)[0]*B(X)[3]-A(X)[2]*B(X)[1])

plt.plot(X,Even)
plt.plot(X,Odd)
plt.plot((0,5),(0,0))
plt.show()
