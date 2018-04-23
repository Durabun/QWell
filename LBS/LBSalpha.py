import cmath
import numpy as np
import math
from scipy.special import *
from numpy.linalg import inv


h = (6.626*10**-34)/(2*np.pi)
m = float(input("mass? "))*9.11*10**-31
E = float(input("energy? "))*9.11*10**-31

n = int(input("number of barriers? "))
S = n+1

V = np.zeros(n+1)
X = np.zeros(n)
beta = np.zeros(n+1)

for l in range(n):
    #The first spot in V corresponds to the free particle, so V0 = 0
    V[l+1] = float(input("Height of barrier "+str(l+1)+": "))*1.602*10**-19
    X[l] = float(input("Thickness of barrier "+str(l+1)+": "))*10**-10
    beta[l+1] = float(input("slope of barrier "+str(l+1)+": "))*10**-9
    l=l+1

#Defining the matrices and K vectors
def K(m,E,V):
    return cmath.sqrt((2*m*(E-V)/h**2))

def Alpha(m,B):
    return ((2*m*B)/h**2)**(1./3)

def Kairy(alpha,x,E,B):
    return cmath.sqrt((alpha**3)*(x-E/B))

def dij(k,E,alpha,B):
    first  = np.matrix(((1,-(1j)/k),(1,(1j)/k)))
    second = np.matrix(((airy(alpha*(-E/B))[0], airy(alpha*(-E/B))[2]),(alpha*airy(alpha*(-E/B))[1],alpha*airy(alpha*(-E/B))[3])))
    return 0.5*np.dot(first,second)

def djk(E,alphaj,xj,Bj,alphak,xk,Bk):
    first = np.matrix(((airy(alphaj*(xj-E/Bj))[0], airy(alphaj*(xj-E/Bj))[2]),(alphakj*airy(alphaj*(xj-E/Bj))[1], alphaj*airy(alphaj*(xj-E/Bj))[3])))

    second = np.matrix(((airy(alphak*(xk-E/Bk))[0], airy(alphak*(xk-E/Bk))[2]),(alphak*airy(alphak*(xk-E/Bk))[1], alphak*airy(alphak*(xk-E/Bk))[3])))

    first = inv(first)
    return np.dot(first, second)

def dkl(k,E,alpha,B):
    second  = np.matrix(((1,-(1j)/k),(1,(1j)*k)))
    first = np.matrix(((airy(alpha*(-E/B))[0], airy(alpha*(-E/B))[2]),(alpha*airy(alpha*(-E/B))[1],alpha*airy(alpha*(-E/B))[3])))
    first = inv(first)
    return np.dot(first,second)

Q = np.matrix(((1,0),(0,1)))
i = 0
j = 1

M = dij(K(m,E,V[0]),E,Alpha(m,beta[0]),beta[0])

for i in range(S):
    if i == n-1:
        j = 0
        N = dkl(K(m,E,V[j]),E,Alpha(m,beta[i]),beta[i])
        break
    
    L = djk(E,Alpha(m,beta[i]),X[i],beta[i],Alpha(m,beta[j]),X[j],beta[j])
    Q = np.dot(Q,M)
    M = np.dot(Q,L)
    i = i+1
    j = j + 1

mat = np.dot(M,N)
print(mat)    

 
