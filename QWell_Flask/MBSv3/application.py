import numpy as np
from flask import Flask, render_template, request
#import jinja2
from MBScalc import *

#Initialize Flask App
app = Flask(__name__)


#First route, renders the number input page
@app.route('/', methods = ['POST','GET'])
def init():
	return render_template("init.html")

#Second route, gets the number, renders the page with remaining inputs
@app.route('/barrier', methods = ['POST','GET'])
def handle():
	number = request.form['number']
	return render_template("barrier.html", number = int(number))

#Final route, handles rest of outputs as well as output
@app.route('/handle', methods = ['POST','GET'])
def barrier():
	#Converts Mass and Energy in SI units
	mass = float(request.form['mass'])*9.11*10**-31
	energy = float(request.form['energy'])*1.602*10**-19
	#Creates temp arrays for barrier heights and thicknesses
	H = request.form.getlist('height')
	D = request.form.getlist('thickness')
	#Loop converts the arrays in floats
	for i in range(len(H)):
		H[i] = float(H[i])
		D[i] = float(D[i])
		i = i + 1
	#Defining the number of boundaries, S
	N = int(len(H))
	S = 2*(N+1)
	#Empty arrays for energy and thickneses
	V = np.zeros(N+1)
	X = np.zeros(N)
	
	#Loop for filling them in
	#These arrays will now be used
	for i in range(N):
		V[i+1] = H[i]*1.602*10**-19
		X[i] = D[i]*10**-10
		i = i +1
	DX = LenBar(X)

	#Calls fn that calculates the tunnel/reflect coefficients
	#Refer to that for more info
	#It returns a tuple
	matrix = CalculateCoeff(mass,energy,N,S,V,X)
	T = matrix[0]
	R = matrix[1]
	
	#Fn that creates a plot file that is displayed
	plot = plots(T,R,K(mass,energy,V[0]),DX)

	#Final render with mass, energy, coeffs, end plot
	return render_template("final.html", mass = mass, energy = energy, out1 = T, out2 = R, image_data = plot)


if __name__=="__main__":
	app.run()
