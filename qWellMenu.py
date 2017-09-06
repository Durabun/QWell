import numpy as np
from os import system
from scipy.optimize import*
import matplotlib.pyplot as plt
import quantumFn
import kivy
#kivy.require('1.0.6') # replace with your current kivy version !

from kivy.app import App
from kivy.uix.label import Label
from kivy.uix.gridlayout import GridLayout
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.button import Button
from kivy.uix.widget import Widget
from kivy.uix.textinput import TextInput
from kivy.uix.anchorlayout import AnchorLayout
from kivy.uix.popup import Popup

h = 6.626*10**(-34)
h_bar = h/(2*np.pi)

class QwellMenu(App):
	
	def build(self):
		
		layout = GridLayout(cols = 2)
		layout.add_widget(Label(text = 'Mass'))
		self.Mass = TextInput(multiline = False)
		layout.add_widget(self.Mass)
		layout.add_widget(Label(text = 'Depth'))
		self.Depth = TextInput(multiline = False)
		layout.add_widget(self.Depth)
		layout.add_widget(Label(text = 'Width'))
		self.hWidth = TextInput(multiline = False)
		layout.add_widget(self.hWidth)
		self.popup = Popup(title = 'Test popup')
		GoBtn = Button(text = 'GO!')
		GoBtn.bind(on_release = self.main)
		layout.add_widget(GoBtn)
		anchor = AnchorLayout(anchor_x = 'center', anchor_y = 'bottom')
		btn = Button(text = ('test'))
		anchor.add_widget(btn)	
		
		return layout

	#This is the main function
	#That consists of most of the code
	def main(self,btn):
		m = (9.11*10**(-31))*float(self.Mass.text)
		V = (1.602*10**(-19))*float(self.Depth.text)
		a = (10**(-10))*float(self.hWidth.text)		
		R = quantumFn.radius(m,a,V)
		x = 0
		i = 1
		j = 0

		while R > i*np.pi:
			i = i+1
		while R > (1+2*j)*np.pi/2:
			j = j+1
		S = i+j

		print (str(S)+" root(s) in total.")

		#Initializes the "x" coord of the even and odd roots
		evenRoot = np.zeros(i)
		oddRoot = np.zeros(j)
		#Used for indexing for the roots
		ie = 0
		jo = 0
		#This is an interative root finder
		#With an x0 used as an initial guess for the Newton function

		dx = 0.0001 #This value has to be smaller to resolve bigger
		x0 = 0.0000  #Scenarios!

		#Used the Even/Odd root functions to establish a function
		#With the radius as the parameter and "x" as a variable
		ERoot = quantumFn.Evenroot(R)
		ORoot = quantumFn.Oddroot(R)

		#Since there will always be one root, this while loop
		#Searches for the first even root, and then puts it into
		#The evenRoot array
		while np.absolute(ERoot(x0)) > 0.05: #This value has to
			x0 = x0 + dx		     #Be small to resolve
		print (x0)			     #Smaller scenarios!
		print (R)
		evenRoot[ie] = newton(ERoot,x0)

		#If there is another even root, it will go into index 1
		#For the even root array
		#The next possible root will be when the radius is at least
		#pi/2, so that is why x0 jumps to that value
		ie = 1
		x0 = (1+2*jo)*np.pi/2

		#This loop will search for the remaining possible roots
		#It goes in the order of looking for an odd root
		#And then an even root.
		while x0 < R:
			while np.absolute(ORoot(x0)) > 0.05:
				x0 = x0 + dx
			oddRoot[jo] = newton(ORoot,x0)
			jo = jo +1
			x0 = ie*np.pi #Resets the initial guess to next possible value
			if x0 < R:
				while np.absolute(ERoot(x0)) > 0.05:
					x0 = x0 +dx
				evenRoot[ie] = newton(ERoot,x0)
				ie = ie +1
				x0 = (1+2*jo)*np.pi/2 #Resets again

		#Initialize the propagation and tunneling vectors
		#At this point I consolidate the even and odd roots
		#Into single arrays instead of having arrays for 
		#Even and odd roots
		propVector = np.zeros(S)
		tunnVector = np.zeros(S)

		ev = 0
		od = 1
		Roots = np.zeros(S)
		#This is to define the vectors of each function
		while ev < S:
			Roots[ev] = evenRoot[int(ev/2)]
			ev = ev +2

		while od < S:
			Roots[od] = oddRoot[int(od/2)]
			od = od + 2

		print (Roots)

		ev = 0
		od = 1

		while ev < S:
			propVector[ev] = (10**(-10))*Roots[ev]/a
			tunnVector[ev] = (10**(-10))*(quantumFn.radFunction(Roots[ev],R))/a
			ev = ev + 2

		while od < S:
			propVector[od] = (10**(-10))*Roots[od]/a
			tunnVector[od] = (10**(-10))*(quantumFn.radFunction(Roots[od],R))/a
			od = od + 2

		print ("Propagation 'k' vectors")
		print (propVector)
		print ("Tunneling 'K' vectors")
		print (tunnVector)

		Bcoeff = np.zeros(S)
		StateEnergy = np.zeros(S)

		l = 0
		while l < S:
			Bcoeff[l] = quantumFn.Beven(a*10**(10),propVector[l],tunnVector[l])
			StateEnergy[l] = -1*quantumFn.energy(quantumFn.radFunction(Roots[l],R),R,V)
			l = l + 2
		l = 1
		while l < S:
			Bcoeff[l] = quantumFn.Bodd(a*10**(10),propVector[l],tunnVector[l])
			StateEnergy[l] = -1*quantumFn.energy(quantumFn.radFunction(Roots[l],R),R,V)
			l = l + 2

		print ("Coefficients")
		print (Bcoeff)
		print ("Energies")
		print (StateEnergy)

		z = 0
		r1 = np.arange((-2*a*10**(10)),(-a*10**(10)),0.0001)
		r2 = np.arange((-a*10**(10)),(a*10**(10)),0.0001)
		r3 = np.arange((a*10**(10)),(2*a*10**(10)),0.0001)
		color = ["r","b","g","c","m","y","k"] #This is the color array

		c = 0
		while z < S:
			if c == 6:
				c = 0
			if z%2 == 1:
				#Note that this file plots the density functions
				#Of the states
				plt1 = quantumFn.op1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
				plt2 = quantumFn.op2(Bcoeff[z],propVector[z])
				plt3 = quantumFn.op3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
		
				plt.plot(r1,plt1(r1)+StateEnergy[z],color[c],r2,plt2(r2)+StateEnergy[z],color[c],r3,plt3(r3)+StateEnergy[z],color[c])
		
		

			else:
				plt1 = quantumFn.ep1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))
				plt2 = quantumFn.ep2(Bcoeff[z],propVector[z])
				plt3 = quantumFn.ep3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**(10)))

				plt.plot(r1,plt1(r1)+StateEnergy[z],color[c],r2,plt2(r2)+StateEnergy[z],color[c],r3,plt3(r3)+StateEnergy[z],color[c])
		
		
			z = z + 1
			c = c + 1

		plt.xlabel("Distance (in Angstrom)")

		return plt.show()


if __name__ == '__main__':
    QwellMenu().run()
