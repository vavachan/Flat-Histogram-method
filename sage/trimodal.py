import numpy as np
import math
from cmath import sqrt
import sys
def g_0(x,p,a):
	return -1/2*(1j*sqrt(3) + 1)*(1/54*((9*p**2*x + 2*x**3 - 9*p*x)*a**2 - 27*p**2*x + 9*p*x)/(a**2*p**3) + 1/6*sqrt(-4/3*a**6*p**4 - 4/3*a**4*x**4 + 4*(a**6 - a**4)*p**3 - 4*(a**6 - 2*a**4 + a**2)*p**2 - 1/3*(a**6 - 2*a**4 + (a**6 + 18*a**4 - 27*a**2)*p**2 + a**2 - 2*(a**6 + 8*a**4 - 9*a**2)*p)*x**2 + 4/3*(a**6 - 3*a**4 + 3*a**2 - 1)*p)/(a**3*p**2))**(1/3) + 1/3*x/p - 1/18*((3*p**2 + x**2 - 3*p)*a**2 + 3*p)*(-1j*sqrt(3) + 1)/(a**2*p**2*(1/54*((9*p**2*x + 2*x**3 - 9*p*x)*a**2 - 27*p**2*x + 9*p*x)/(a**2*p**3) + 1/6*sqrt(-4/3*a**6*p**4 - 4/3*a**4*x**4 + 4*(a**6 - a**4)*p**3 - 4*(a**6 - 2*a**4 + a**2)*p**2 - 1/3*(a**6 - 2*a**4 + (a**6 + 18*a**4 - 27*a**2)*p**2 + a**2 - 2*(a**6 + 8*a**4 - 9*a**2)*p)*x**2 + 4/3*(a**6 - 3*a**4 + 3*a**2 - 1)*p)/(a**3*p**2))**(1/3))
def taylor(x,p,a):
	return ((243*1j*sqrt(3)*(3*(-1)**(1/3) - 14) - 243*1j*sqrt(3)*(3*(-1)**(1/3) - 20) - 1458*1j*sqrt(3))*a**8 + (-729*1j*sqrt(3)*(3*(-1)**(1/3) - 14) + 729*1j*sqrt(3)*(3*(-1)**(1/3) - 20) + 4374*1j*sqrt(3))*a**6 + (729*1j*sqrt(3)*(3*(-1)**(1/3) - 14) - 729*1j*sqrt(3)*(3*(-1)**(1/3) - 20) - 4374*1j*sqrt(3))*a**4 + ((-243*1j*sqrt(3)*(3*(-1)**(1/3) - 14) + 243*1j*sqrt(3)*(3*(-1)**(1/3) - 20) + 1458*1j*sqrt(3))*a**8 + (6561*1j*sqrt(3)*(5*(-1)**(1/3) + 14) - 59049*1j*sqrt(3)*((-1)**(1/3) + 2) - 26244*(-1)**(1/3) + 26244)*a**4 + (-6561*1j*sqrt(3)*(5*(-1)**(1/3) + 14) + 59049*1j*sqrt(3)*((-1)**(1/3) + 2) + 26244*(-1)**(1/3) - 26244)*a**2)*p**3 + (-243*1j*sqrt(3)*(3*(-1)**(1/3) - 14) + 243*1j*sqrt(3)*(3*(-1)**(1/3) - 20) + 1458*1j*sqrt(3))*a**2 + ((729*1j*sqrt(3)*(3*(-1)**(1/3) - 14) - 729*1j*sqrt(3)*(3*(-1)**(1/3) - 20) - 4374*1j*sqrt(3))*a**8 + (729*1j*sqrt(3)*(27*(-1)**(1/3) + 98) - 729*1j*sqrt(3)*(27*(-1)**(1/3) + 44) - 39366*1j*sqrt(3))*a**6 + (2187*1j*sqrt(3)*(37*(-1)**(1/3) + 8) - 2187*1j*sqrt(3)*(25*(-1)**(1/3) - 4) + 26244*(-1)**(1/3) - 26244)*a**4 + (6561*1j*sqrt(3)*(5*(-1)**(1/3) + 14) - 59049*1j*sqrt(3)*((-1)**(1/3) + 2) - 26244*(-1)**(1/3) + 26244)*a**2)*p**2 + ((-729*1j*sqrt(3)*(3*(-1)**(1/3) - 14) + 729*1j*sqrt(3)*(3*(-1)**(1/3) - 20) + 4374*1j*sqrt(3))*a**8 + (-729*1j*sqrt(3)*(9*(-1)**(1/3) + 70) + 729*1j*sqrt(3)*(9*(-1)**(1/3) + 52) + 13122*1j*sqrt(3))*a**6 + (729*1j*sqrt(3)*(27*(-1)**(1/3) + 98) - 729*1j*sqrt(3)*(27*(-1)**(1/3) + 44) - 39366*1j*sqrt(3))*a**4)*p)*x**3/(104976*1j*(-1)**(1/6)*a**8*p**6 + (-419904*1j*(-1)**(1/6)*a**8 + 419904*1j*(-1)**(1/6)*a**6)*p**5 + (629856*1j*(-1)**(1/6)*a**8 - 1259712*1j*(-1)**(1/6)*a**6 + 629856*1j*(-1)**(1/6)*a**4)*p**4 + (-419904*1j*(-1)**(1/6)*a**8 + 1259712*1j*(-1)**(1/6)*a**6 - 1259712*1j*(-1)**(1/6)*a**4 + 419904*1j*(-1)**(1/6)*a**2)*p**3 + (104976*1j*(-1)**(1/6)*a**8 - 419904*1j*(-1)**(1/6)*a**6 + 629856*1j*(-1)**(1/6)*a**4 - 419904*1j*(-1)**(1/6)*a**2 + 104976*1j*(-1)**(1/6))*p**2) - 1/3888*((27*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*(2*(-1)**(1/3) - 5) + 27*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) + 8) + 243*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a**5 - (54*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*(2*(-1)**(1/3) - 5) + 54*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) + 8) + 486*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a**3 + ((27*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*(2*(-1)**(1/3) - 5) + 27*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) + 8) + 243*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a**5 + (486*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1) - 324*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 2) + 486*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a**3 - (729*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1) - 486*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 2) + 729*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a)*p**2 + (27*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*(2*(-1)**(1/3) - 5) + 27*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) + 8) + 243*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a - ((54*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*(2*(-1)**(1/3) - 5) + 54*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) + 8) + 486*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a**5 - (486*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1) - 324*sqrt(3)*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 2) + 486*1j*(a**2*p - a**2 + 1)**(2/3)*p**(2/3)*((-1)**(1/3) - 1))*a)*p)*x**2/((-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**6*p**(31/6) - 3*((-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**6*p**(1/6) - (-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**4*p**(1/6))*p**4 + 3*((-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**6*p**(1/6) - 2*(-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**4*p**(1/6) + (-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**2*p**(1/6))*p**3 - ((-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**6*p**(1/6) - 3*(-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**4*p**(1/6) + 3*(-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*a**2*p**(1/6) - (-1)**(1/6)*(a**2*p - a**2 + 1)**(1/6)*p**(1/6))*p**2) - ((-9*1j*sqrt(3)*((-1)**(1/3) + 1) - 9*(-1)**(1/3) + 36*1j*(-1)**(1/6) + 9)*a**2 + ((9*1j*sqrt(3)*((-1)**(1/3) + 1) + 9*(-1)**(1/3) - 36*1j*(-1)**(1/6) - 9)*a**2 - 27*1j*sqrt(3)*((-1)**(1/3) + 1) - 27*(-1)**(1/3) + 27)*p + 9*1j*sqrt(3)*((-1)**(1/3) + 1) + 9*(-1)**(1/3) - 36*1j*(-1)**(1/6) - 9)*x/(108*1j*(-1)**(1/6)*a**2*p**2 + (-108*1j*(-1)**(1/6)*a**2 + 108*1j*(-1)**(1/6))*p) + 1/18*sqrt(3)*(-3*1j*sqrt(3)*sqrt(a**2*p - a**2 + 1)*((-1)**(5/6)*p**(1/6) + (-1)**(1/6)*p**(1/6)) + 3*sqrt(a**2*p - a**2 + 1)*((-1)**(5/6)*p**(1/6) - (-1)**(1/6)*p**(1/6)))/(a*p**(2/3))	
def bimodal_1(x,a):
	return 1/2*(a**2 - sqrt(a**4 + 4*a**2*x**2 - 2*a**2 + 1) - 1)/(a**2*x)
def bimodal_2(x,a):
	return 1/2*(a**2 + sqrt(a**4 + 4*a**2*x**2 - 2*a**2 + 1) - 1)/(a**2*x) 
def bimodal_3(x,a):
	return np.log(-(a*x + sqrt((a**2 - 1)*x**2 + 1))/(x - 1))
def bimodal_4(x,a):
	return np.log(-(a*x - sqrt((a**2 -1)*x**2 + 1))/(x - 1))
def f(x,p,beta,h):
	return p*np.log(np.cosh(x))+(1-p)/2*np.log(np.cosh(x+beta*h)*np.cosh(x-beta*h)) 
def F(x,beta):
	return beta*x**2/2
def f_star(x,p,beta,h):
	x_star=np.arctanh(g_0(x,p,np.tanh(beta*h)).real)
	#print x,g_0(x,p,np.tanh(beta*h)).real
	return x*x_star-f(x_star,p,beta,h) 
def I(x,p,beta,h):
	return f_star(x,p,beta,h)-F(x,beta)
def der_g(x,p,beta,h):
	return (np.arctanh(g_0(x+.005,p,np.tanh(beta*h)).real)-np.arctanh(g_0(x,p,np.tanh(beta*h)).real))/.005
p=float(sys.argv[1])
h=float(sys.argv[2])
beta=1.0/float(sys.argv[3])
def a(beta,h):
	return np.tanh(beta*h)
b=a(beta,h)
#c=np.cosh(beta*h)
min=I(-0.999,p,beta,h)
x_min=-1
for i in range (1,2000-1):
	x=(i-1000)/1000.0
	if(I(x,p,beta,h)<min):
		min=I(x,p,beta,h)	
		x_min=x
########print(abs(x_min),succ(x_min,p,beta,h))
#print(1.0/beta,abs(x_min),beta/(-beta+der_g(x_min,p,beta,h)))
print(abs(x_min))
#for i in range (1,2000-1):
	#x=(i-1000)/1000.0
	#print(x,f_star(x,p,beta,h)-F(x,beta)-min,g_0(x,p,b).real,np.tanh(beta*x),abs(x_min))#I(0,p,beta,h)-min)#f(x,p,beta,h),x*g_0(x,p,b).real)#np.tanh(beta*x),g_0(x,p,b).real,g_1(x,p,b).real,g_2(x,p,b).real)