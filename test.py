'''
This is a script to test out the fsolve function. It seems
to work pretty well so I may end up using it.
The idea is to subtract 2 functions, then use fsolve
to find the root of the difference. That root will be were the two functions intersect.
'''


import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

#x = np.linspace(0,5,1000)
#y1 = x*(np.sin(x))/(np.cos(x))
#y2 = -x*(np.cos(x))/(np.sin(x))
#yr = np.sqrt(10-x**2)

def diff1(u):
	return  u*(np.sin(u))/(np.cos(u))- np.sqrt(49-u**2)
def diff2(u):
	return -u*(np.cos(u))/(np.sin(u))- np.sqrt(49-u**2)


root =fsolve(diff1,1)
root2 =fsolve(diff1,3.1)
root3 = fsolve(diff1,5)


print root
print root2
print root3
print fsolve(diff2,2)
print fsolve(diff2,4)

#plt.ylim(0,10)
#plt.plot(x,y1)
#plt.plot(x,y2)
#plt.plot(x,yr)
#plt.grid()
#plt.show()
