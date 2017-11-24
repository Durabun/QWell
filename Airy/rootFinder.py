import numpy as np
import matplotlib.pyplot as plt
from scipy.special import*
from scipy.optimize import*


r = np.arange(-10,0,0.1)
X = np.arange(0.1,4,0.1)
x = 0
S = 0


V = 1.602*10**(-19)
a = 0.5*10**(-10)
beta = V/a
m = 206*9.11*10**(-31)
hb = (6.626*10**(-34))/(2*np.pi)
alpha = ((2*m*beta)/hb**2)**(1./3)


def Even(alpha,a,x):
	A = airy(-x)
	B = airy(alpha*a-x)
	K = (a*alpha**3-alpha**2*x)**(1./2)

	return K*(A[1]*B[2]-A[3]*B[0])+alpha*(A[1]*B[3]-A[3]*B[1])

def Odd(alpha,a,x):
	A = airy(-x)
	B = airy(alpha*a-x)
	K = (a*alpha**3-alpha**2*x)**(1./2)

	return K*(A[0]*B[2]-A[3]*B[0])+alpha*(A[0]*B[3]-A[2]*B[1])

x0 = 0.5
dx = 0.01

while(np.sign(Even(alpha,a,x0)) == np.sign(Even(alpha,a,x0+dx))):
	print(x0)
	x0 = x0+dx

print(Even(alpha,a,x0))
#print(Even(alpha,a,x0))

plt.plot(X,Even(alpha,a,X))
plt.plot(X,Odd(alpha,a,X))
#plt.ylim(-1,1)

plt.show()
