import numpy as np
from os import system
from scipy.optimize import*

#For a one dimensional well...

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
	return np.absolute(V)*(y/r)**2
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
	I	|  	 II   |	  III
		|         |
		|    	  |
	  -a|____0____|a

'''
	
	
	
	
#The following set of functions define the wavefunction
#In each of the three regions, with an "e" or "o"
#Denoting an even or odd function for the respective root
def ef1(B,K,k,a):
	return lambda x: 2*B*np.exp((a+x)*K)*np.cos(a*k)
	
def ef2(B,k):
	return lambda x: 2*B*np.cos(k*x)	
	
def ef3(B,K,k,a):
	return lambda x: 2*B*np.exp((a-x)*K)*np.cos(a*k)	
	
def of1(B,K,k,a):
	return lambda x: -2*B*np.exp((a+x)*K)*np.sin(a*k)
	
def of2(B,k):
	return lambda x: 2*B*np.sin(k*x)	
	
def of3(B,K,k,a):
	return lambda x: 2*B*np.exp((a-x)*K)*np.sin(a*k)