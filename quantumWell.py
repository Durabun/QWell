import numpy as np
from os import system
from scipy.optimize import*
import matplotlib.pyplot as plt
import quantumFn

h = 6.626*10**(-34) #These are some constants
h_bar = h/(2*np.pi) #Specifically, Planck's Constant

m = float(input("Mass of particle: "))
V = float(input("Depth of well (J): "))
a = float(input("Half width of well centered at x=0: "))

#R = quantumFn.radius(m,a,V)
R = m+V+a
print (R)