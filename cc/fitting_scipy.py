import numpy as np
import scipy
from scipy import optimize
import math
import sys
m=float(sys.argv[1])
def func(x,a,b,c):
	return a*x**2+b*x+c
exec "M=np.loadtxt('m_%d.dat')" %m
A=M[:,2]
B=M[:,0]
popt,pcov = optimize.curve_fit(func,B,A)
print len(A)
for i in range(0,len(A)):
	print B[i],A[i],func(B[i],popt[0],popt[1],popt[2])
