import numpy as np
from scipy.optimize import *
from flask import Flask, render_template, request
from quantumFn import *
import matplotlib.pyplot as plt
import os, time, glob

#This will port the Kivy app to a flask app

h_bar = (6.626*10**(-34))/(2*np.pi)

app = Flask(__name__)

@app.route('/', methods = ['POST','GET'])
def init():
    return render_template("sqInput.html")

@app.route('/square', methods = ['POST','GET'])
def sqCalc():
    m = float(request.form["mass"])*(9.11*10**(-31))
    V = float(request.form["depth"])*(1.602*10**(-19))
    a = float(request.form["width"])*(10**(-10))

    R = radius(m,a,V)
    i = 1
    j = 0

    while R > i*np.pi:
        i = 1+1
    while R > (1+2*j)*np.pi/2:
        j = j+1
    S = i+j

#    print(str(S)+" root(s) in total.")

    evenRoot = np.zeros(i)
    oddRoot = np.zeros(j)

    ie = 0
    jo = 0

    dx = 0.0001

    x0 = 0.0001

    ERoot = Evenroot(R)
    ORoot = Oddroot(R)

    while np.absolute(ERoot(x0)) > 0.05:
        x0 = x0 + dx

#    print(x0)
#    print(R)
    evenRoot[ie] = newton(ERoot, x0)	
    
    ie = 1
    x0 = (1+2*jo)*np.pi/2

    while x0 < R:
        while np.absolute(ORoot(x0)) > 0.05:
            x0 = x0 + dx
        oddRoot[jo] = newton(ORoot,x0)
        jo = jo +1
	x0 = ie*np.pi

	if x0 < R:
            while np.absolute(ERoot(x0)) > 0.05:
		x0 = x0 +dx
	    evenRoot[ie] = newton(ERoot,x0)
	    ie = ie+1
	    x0 = (1+2*jo)*np.pi/2

    propVector = np.zeros(S)
    tunnVector = np.zeros(S)
    Bcoeff = np.zeros(S)
    StateEnergy = np.zeros(S)

    ev = 0
    od = 1
    Roots = np.zeros(S)

    while ev < S:
        Roots[ev] = evenRoot[int(ev/2)]
	ev = ev+2

    while od < S:
	Roots[od] = oddRoot[int(od/2)]
	od = od+2

    ev = 0
    od = 1

    while ev < S:
        propVector[ev] = (10**-10)*Roots[ev]/a
	tunnVector[ev] = (10**-10)*(radFunction(Roots[ev],R))/a
	ev = ev+2

    while od < S:
	propVector[od] = (10**-10)*Roots[od]/a
	tunnVector[od] = (10**-10)*(radFunction(Roots[od],R))/a
	od = od+2

    l = 0

    while l < S:
        Bcoeff[l] = Beven(a*10**10,propVector[l],tunnVector[l])
	StateEnergy[l] = -1*energy(radFunction(Roots[l],R),R,V)
        l = l +2

    l = 1
    while l < S:
        Bcoeff[l] = Bodd(a*10**10,propVector[l],tunnVector[l])
	StateEnergy[l] = -1*energy(radFunction(Roots[l],R),R,V)
	l = l+1

    plot = Plots(a,S,Bcoeff,tunnVector,propVector,StateEnergy)
#    z = 0
#    r1 = np.arange((-2*a*10**10),(-a*10**10),0.0001)
#    r2 = np.arange((-a*10**10),(a*10**10),0.0001)
#    r3 = np.arange((a*10**10),(2*a*10**10),0.0001)
#    color = ['r','b','g','c','m','y','k']

#    plt.figure()
#    axes = plt.gca()
#    plt.title("Finite Square Well")
#    plt.xlabel("Angstrom")
#    plt.ylabel("Probability Amplitude")
#    c = 0
#    while z < S:
#        if c == 6:
#            c = 0
#        if z%2 == 1:
#            plt1 = op1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**10))
#            plt2 = op2(Bcoeff[z],propVector[z])
#            plt3 = op3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**10))
	    
#            plt.plot(r1,plt1(r1)+StateEnergy[z],color[c],r2,plt2(r2)+StateEnergy[z],color[c],r3,plt3(r3)+StateEnergy[z],color[c])


#        else:
#            plt1 = ep1(Bcoeff[z],tunnVector[z],propVector[z],(a*10**10))
#            plt2 = ep2(Bcoeff[z],propVector[z])
#            plt3 = ep3(Bcoeff[z],tunnVector[z],propVector[z],(a*10**10))
	    
#            plt.plot(r1,plt1(r1)+StateEnergy[z],color[c],r2,plt2(r2)+StateEnergy[z],color[c],r3,plt3(r3)+StateEnergy[z],color[c])


#        z = z +1
#	c = c +1


#    if not os.path.isdir('static'):
#        os.mkdir('static')
#    else:
#        for filename in glob.glob(os.path.join('static','*.png')):
#            os.remove(filename)
#    plotfile = os.path.join('static', str(time.time())+'.png')
#    plt.savefig(plotfile)

#    plt.show()
	
    return render_template("sqDisplay.html", S = S, Energy = StateEnergy, Roots = Roots, PVector = propVector, TVector = tunnVector, PLOT = plot)









if __name__=="__main__":
	app.run()
