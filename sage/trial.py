import numpy as np
import sys
p=float(sys.argv[1])
h=float(sys.argv[2])
def A(p,a):
	return p+(1-a**2)*(1-p)
for t in range (50,300):
	T=t*1.0/100
	beta=1.0/T
	a=np.tanh(beta*h)
	print T,beta*a**2*A(p,a)-p*a**3-beta**3*A(p,a)**3/3,p*a*beta**3-a**2*beta**2+A(p,a)/3
