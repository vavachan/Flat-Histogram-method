import numpy as np
import scipy
from scipy import optimize
import math
import sys
L=float(sys.argv[1])
def func(x,p,a,b,c,d):
	#print x
	return p*x**4+b*x**2+d
exec "A=np.loadtxt('data_%d.dat')" %L
N=int(L*L*L)
entropy_L=np.zeros(shape=(N+1,3*N+1))
entropy=np.zeros(shape=(N+1,3*N+1))
for m in range (0,N+1):
        for n in range (int(A[m][2]),int(A[m][1]),2):
                entropy_L[m][n]=1
#now entropy_L2 has non-zero value if the configuration class is allowed.
entropy_L[0][3*N]=1
entropy_L[1][3*N-6]=1
entropy_L[N][3*N]=1
entropy_L[N-1][3*N-6]=1
exec "A=np.loadtxt('entropy_mp_9.dat')" 
for i in range (0,len(A)):
        entropy[int((A[i][1]+N)/2)][int(A[i][0])]=float(A[i][2])


for n in range (400,402,2):
	A=[]
	B=[] #n=int(sys.argv[2])
	for m in range (0,N+1):
		if(entropy[m][n]!=0):
			#print entropy[m][n],m,n
			A.append((entropy[m][n]))
			B.append(m)
	#print A,B
	X=np.zeros(shape=len(A))
	Y=np.zeros(shape=len(A))
	for i in range (len(A)):
		X[i]=2*B[i]-N
		Y[i]=A[i]
	popt,pcov = optimize.curve_fit(func,X,Y)
	for m in range (0,N+1):
		if(entropy_L[m][n]!=0 ):
		#	break
			print n,m,func(2*m-N,popt[0],popt[1],popt[2],popt[3],popt[4]),entropy[m][n]
			entropy[m][n]=func(2*m-N,popt[0],popt[1],popt[2],popt[3],popt[4])
for n in range (0,3*N+1):
	for m in range (0,N+1):
		if(entropy[m][n]!=0):
			break
			print 2*m-N,n,entropy[m][n]
