import numpy as np
from flask import Flask, render_template, request
import jinja2
from MBScalc import *

#Init Flask App
app = Flask(__name__)

@app.route('/', methods = ['GET','POST'])
def main():
	number  = int(request.form['number'])
	if request.method == 'POST':
		result = number
	else:
		result = 0
#	mass = float(request.form['mass'])*9.11*10**-31
#	energy = float(request.form['energy'])*1.602*10**-19
#
#	H = request.form.getlist('height')
#	D = request.form.getlist('thickness')
#
#	for i in range(len(H)):
#		H[i] = float(H[i])
#		D[i] = float(D[i])
#		i = i +1
#
#	N = int(len(H))
#	S = 2*(N+1)
#
#	V = np.zeros(N+1)
#	X = np.zeros(N)

#	for i in range(N):
#		V[i+1] = H[i]*1.602*10**-19
#		X[i] = D[i]*10**-10
#		i = i + 1
#	DX = LenBar(X)
#
#	matrix = Calculateoeff(mass,energy,N,S,V,X)
#	T = matrix[0]
#	R = matrix[1]
#
#	plot = plots(T,R,K(mass,energy,V[0]),DX)
#
	return render_template("view.html",num = number)#, mass = mass, energy = energy, number = number, out1 = T, out2 = R, image_data = plot)

if __name__=="__main__":
	app.run()
