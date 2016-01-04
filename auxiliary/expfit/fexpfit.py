#!/usr/bin/env python
# I read s [slab thickness], a [<sig^2>/<sig^2>], and lamcsig [log-normal lamc],
# generate the real log-normal covariance function, and curve fit an exponential
# to that data, yielding lamcw [Gaussian lamcw].  This allows use of exponential
# covariance with fitted correlation length with only minimal error to the
# real covariance.
import numpy as np
from scipy.optimize import curve_fit

def origfunc(x,lamcsig,a):  #func we are fitting other to
    return np.log(np.exp(-x/lamcsig)*(a-1.0)+1.0)/np.log(a)

def fittingfunc(x,lamcw):  #func we are using to fit for param lamcw
    return np.exp(-x/lamcw)

#read input
s,a,lamcsig = np.loadtxt('auxiliary/expfit/s.a.lamcsig.txt')

x = np.linspace(0.0000001,0.0000001+s,100)
yn= origfunc(x,lamcsig,a)
lamcw,pcov = curve_fit(fittingfunc,x,yn)  #fit curve

#print output
np.savetxt('auxiliary/expfit/lamcw.txt',lamcw)
