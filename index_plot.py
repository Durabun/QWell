import numpy as np
import matplotlib.pyplot as plt

evenR = np.array([1.212,3.368])
oddR = np.array([2.381])
S = evenR.size+oddR.size
prop = np.zeros(S)
tunn = np.zeros(S)
i=0
j=1
a=0.469

def rad(x):
	return np.sqrt((1.1*np.pi)**2-x**2)
	
print (S)

print (prop)
print (tunn)

while i< evenR.size:
	prop[2*i] = evenR[i]/a
	tunn[2*i] = rad(evenR[i])/a
	i=i+1
	print (i)
print("odd")
while j-1 < oddR.size:
	prop[j] = oddR[j-1]/a
	tunn[j] = rad(oddR[j-1])/a
	j=j+2
	print (j)

print (prop)
print (tunn)

Bcoeff = np.array([0.6318,0.6171,0.4823])

#def Bfn(k,K):
#	return k+1.2*K
	
l = 0

#while l < S:
#	Bcoeff[l] = Bfn(prop[l],tunn[l])
#	l=l+1

print (Bcoeff)

z = 0
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

	
r1 = np.arange(-5,-a,0.001)
r2 = np.arange(-a,a,0.001)
r3 = np.arange(a,5,0.001)
color = ["r","b","g"]
while z <S:
#	plt.figure
	if z%2 == 1:
		plt1 = of1(Bcoeff[z],tunn[z],prop[z],a)
		plt2 = of2(Bcoeff[z],prop[z])
		plt3 = of3(Bcoeff[z],tunn[z],prop[z],a)
		plt.plot(r1,plt1(r1),color[z],r2,plt2(r2),color[z],r3,plt3(r3),color[z])
#		plt.plot(r2,plt2(r2))
#		plt.plot(r3,plt3(r3))
	else:
		plt1 = ef1(Bcoeff[z],tunn[z],prop[z],a)
		plt2 = ef2(Bcoeff[z],prop[z])
		plt3 = ef3(Bcoeff[z],tunn[z],prop[z],a)
		plt.plot(r1,plt1(r1),color[z],r2,plt2(r2),color[z],r3,plt3(r3),color[z])
#		plt.plot(r2,plt2(r2))
#		plt.plot(r3,plt3(r3))
	z = z+1
		
plt.show()