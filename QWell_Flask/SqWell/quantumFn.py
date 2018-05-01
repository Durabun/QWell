#QWell 0.5, the quantum well simulator!
#Made by Jonathan Eugenio
#And Mark Romero

#Please refer to the license in the github repo (https://github.com/Durabun/QWell.git)

import numpy as np
from os import system
import matplotlib.pyplot as plt
from scipy.optimize import*
import os, time, glob

#For a one dimensional well...


h = 6.626*10**(-34) #These are some constants
h_bar = h/(2*np.pi) #Specifically, Planck's Constant

#The fn below defines the even root function
def Evenroot(r):
	return lambda x:x*((np.sin(x))/(np.cos(x)))-np.sqrt(r**2-x**2)


#The fn below defines the odd root function
def Oddroot(r):
	return lambda x:-x*((np.cos(x))/(np.sin(x)))-np.sqrt(r**2-x**2)

#Using a lambda function enables x as a variable while changing other variables

#The function below determines the radius of possible roots 
#Based on the parameters of the particle and the potential well
def radius(m,a,V):
	return np.sqrt((2*m*np.absolute(V)*a**2)/h_bar**2)
	
#This function is used to determine the energy of the state
def energy(y,r,V):
	return (1/(1.602*10**(-19)))*np.absolute(V)*(y/r)**2

#This function is used to determine the "Y" value of the root
#Needed for the energy and for the tunneling vector
def radFunction(x,r):
	return np.sqrt(r**2-x**2)

#Defining the propagation "k" vector as a function of the root
def propVector(x,a):
	return (x/a)*10**(-10)

#Defining the propagation "kappa" vector as a function of the
#"Y" value of the root
def tunnelVector(y,a):
	return (y/a)*10**(-10)


#In a quantum well where "a" is the half width of the well
#There are three regions:
#Region I where -infinity <=x<=-a
#Region II where -a<=x<=a
#Region III where a<=x<= infinity
#Because it is a finite well, the wave function exists in each of the three regions
#As opposed to an infinite well, there the function cannot exist in regions I and III
#Because it would not satisfy the time-independent Schrodinger Equation



'''
_________         __________
	|   	  |
  I	|    II   |    III
	|         |
	|    	  |
      -a|____0____|a
'''
	

#B coefficient for the even waves
def Beven(a,k,K):
	return 1/np.sqrt((4*(np.cos(a*k))**2/K)+(2/k)*(2*a*k+np.sin(2*a*k)))
#B coefficient for the odd waves
def Bodd(a,k,K):
	return 1/np.sqrt((4*(np.sin(a*k))**2/K)+4*(a-(np.sin(2*a*k))/(2*k)))
	
	
#The following set of functions define the wavefunction
#In each of the three regions, with an "e" or "o"
#Denoting an even or odd function for the respective root
def ef1(B,K,k,a):
	return lambda x: 2*B*np.exp((a+x)*K)*np.cos(a*k)
def ep1(B,K,k,a):
	return lambda x: (2*B*np.exp((a+x)*K)*np.cos(a*k))**2
	
def ef2(B,k):
	return lambda x: 2*B*np.cos(k*x)	
def ep2(B,k):
	return lambda x: (2*B*np.cos(k*x))**2	
	
def ef3(B,K,k,a):
	return lambda x: 2*B*np.exp((a-x)*K)*np.cos(a*k)	
def ep3(B,K,k,a):
	return lambda x: (2*B*np.exp((a-x)*K)*np.cos(a*k))**2
	
def of1(B,K,k,a):
	return lambda x: -2*B*np.exp((a+x)*K)*np.sin(a*k)

def op1(B,K,k,a):
	return lambda x: (-2*B*np.exp((a+x)*K)*np.sin(a*k))**2
	
def of2(B,k):
	return lambda x: 2*B*np.sin(k*x)	
def op2(B,k):
	return lambda x: (2*B*np.sin(k*x))**2	
	
def of3(B,K,k,a):
	return lambda x: 2*B*np.exp((a-x)*K)*np.sin(a*k)

def op3(B,K,k,a):
	return lambda x: (2*B*np.exp((a-x)*K)*np.sin(a*k))**2

#Flask Plot here
def Plots(a,S, Bcoeff, TunVec, PropVec,Energy):
    z = 0
    r1 = np.arange((-2*a*10**10),(-a*10**10),0.0001)
    r2 = np.arange((-a*10**10),(a*10**10),0.0001)
    r3 = np.arange((a*10**10),(2*a*10**10),0.0001)
    color = ['r','b','g','c','m','y','k']


    plt.figure()
    axes = plt.gca()
    plt.title("Finite Square Well")
    plt.xlabel("Angstrom")
    plt.ylabel("Probability Amplitude")
    c = 0

    while z < S:
        if c == 6:
            c = 0
        if z%2 == 1:
            plt1 = op1(Bcoeff[z],TunVec[z],PropVec[z],(a*10**10))
            plt2 = op2(Bcoeff[z],PropVec[z])
            plt3 = op3(Bcoeff[z],TunVec[z],PropVec[z],(a*10**10))
	    
            plt.plot(r1,plt1(r1)+Energy[z],color[c],r2,plt2(r2)+Energy[z],color[c],r3,plt3(r3)+Energy[z],color[c])


        else:
            plt1 = ep1(Bcoeff[z],TunVec[z],PropVec[z],(a*10**10))
            plt2 = ep2(Bcoeff[z],PropVec[z])
            plt3 = ep3(Bcoeff[z],TunVec[z],PropVec[z],(a*10**10))
	    
            plt.plot(r1,plt1(r1)+Energy[z],color[c],r2,plt2(r2)+Energy[z],color[c],r3,plt3(r3)+Energy[z],color[c])


        z = z +1
	c = c +1

    if not os.path.isdir('static'):
        os.mkdir('static')
    else:
        for filename in glob.glob(os.path.join('static','*.png')):
            os.remove(filename)
    plotfile = os.path.join('static', str(time.time())+'.png')
    plt.savefig(plotfile)
    return plotfile


    
