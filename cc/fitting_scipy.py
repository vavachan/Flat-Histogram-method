import numpy as np
import scipy
from scipy import optimize
import math
import sys
from scipy import interpolate
import pylab as pl
def func(x,a,b,c):
	return a*x**2+b*x+c
exec "M=np.loadtxt('entropy_1para_6.dat')" 
N=6**3
A=M[:,1]/N
B=M[:,0]/N
tck = interpolate.splrep(B,A,s=0)
xnew=np.arange(0,3*8**3,2)
xnew=xnew/(8.0**3)
ynew = interpolate.splev(xnew, tck, der=0)
pl.plot(xnew,ynew)
pl.plot(B,A)
pl.show()
np.savetxt('entropy_1para_ini.dat', np.c_[xnew*(8**3),ynew])

