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

Even = K(X)*(A(X)[1]*B(X)[2]-A(X)[3]*B(X)[0])+alpha*(A(X)[1]*B(X)[3]-A(X)[3]*B(X)[1])

Odd = K(X)*(A(X)[0]*B(X)[2]-A(X)[2]*B(X)[0])+alpha*(A(X)[0]*B(X)[3]-A(X)[2]*B(X)[1])

#Figure out using fsolver/newton's method to find roots!

plt.plot(X,Even)
plt.plot(X,Odd)
plt.show()
