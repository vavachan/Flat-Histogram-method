import numpy as np
import sys
def func(p,h,beta):
#	return np.cosh(beta*h)**2
	return 1.0/(p+(1.0/(np.cosh(beta*h))**2)*(1-p))
def func2(a):
	return 1/3**(0.5)*a/(1-a**2)
def der_f(p,h,beta):
	return (func(p,h,beta+.005)-func(p,h,beta))/.005
def doub_der(p,h,beta):
	return (der_f(p,h,beta+.005)-der_f(p,h,beta))/.005
beta=0.0001	
p=float(sys.argv[1])
h=float(sys.argv[2])
########print (p+(1-(np.tanh(beta*h))**2)*(1-p))
########print doub_der(p,h,beta)
while 1:
	a=np.tanh(beta*h)
	beta_new=func2(a)#(p,h,beta)
	print beta
	if(beta-.00001<beta_new<beta+.00001):
		break
	else:
		beta=beta_new
print 1/beta
