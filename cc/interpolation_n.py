import numpy as np
import scipy
from scipy import optimize
import math
import sys
from scipy import interpolate
L=float(sys.argv[1])
L1=L-2
def func(x,p,a,b,c,d):
	#print x
	return p*x**4+b*x**2+d
def func2(x,a,b,c):
	return a*x**2+b*x+c
exec "A=np.loadtxt('%d.dat')" %L
N=int(L*L*L)
N1=int(L1*L1*L1)
entropy_L=np.zeros(shape=(N+1,3*N+1))
entropy=np.zeros(shape=(N1+1,3*N1+1))
for m in range (0,N+1):
        for n in range (int(A[m][3]),int(A[m][1])+2,2):
                entropy_L[m][n]=1
#now entropy_L2 has non-zero value if the configuration class is allowed.
entropy_L[0][3*N]=1
entropy_L[1][3*N-6]=1
entropy_L[N][3*N]=1
entropy_L[N-1][3*N-6]=1
exec "A=np.loadtxt('entropy_mp_9.dat')" 
for i in range (0,len(A)):
        entropy[int((A[i][1]+N1)/2)][int(A[i][0])]=float(A[i][2])
for n in range (40,540,2):
	A=[]
	B=[] #n=int(sys.argv[2])
	for m in range (0,N1+1):
		if(entropy[m][n]!=0):
			#print entropy[m][n],m,n
			A.append(entropy[m][n]/N1)
			B.append(m*1.0/N1)
	#print A,B
	X=np.zeros(shape=len(A))
	Y=np.zeros(shape=len(A))
	for i in range (len(A)):
		X[i]=2*B[i]-1
		Y[i]=A[i]
	popt,pcov = optimize.curve_fit(func,X,Y)
	n1=int(N*n/N1)
	for m in range (0,N+1):
		if(entropy_L[m][n1]!=0):
			m1=(2.0*m-N)/N
			m2=int(N1*m1)
			#m2=(int(m1)+N1)/2#(N1/N)*(2*m-N)
		#	break
		#	print n,m1,m2*1.0/N1,func(m1,popt[0],popt[1],popt[2],popt[3],popt[4]),entropy[int((m2+N1)/2.0)][n]
			entropy_L[m][n1]=N*func(m1,popt[0],popt[1],popt[2],popt[3],popt[4])
#print int(N*560.0/N1)
for m in range (50,N-50):
	A=[]
	B=[]
	m1=m*1.0/N
	m1=int(m1*N1)
	for n in range (0,3*N+1):
		if (entropy_L[m][n]!=1 and entropy_L[m][n]!=0):
			A.append(entropy_L[m][n])
			B.append(n)		
	xnew=[]
	X=np.zeros(shape=len(A))
	Y=np.zeros(shape=len(A))
	for i in range (len(A)):
		X[i]=B[i]*1.0/N
		Y[i]=A[i]*1.0/N
	tck = interpolate.splrep(X,Y,s=0)
	#print m,int(B[0])
	for n in range (int(B[0])-2,int(B[len(A)-1])+1,2):
		if (entropy_L[m][n]!=0):
			xnew.append(n*1.0/N)		
	#tck = interpolate.splrep(y,x,s=0)				
	ynew = interpolate.splev(xnew, tck, der=0)	
	for n in range (0,len(ynew)):
	#	break
		entropy_L[m][int(N*xnew[n])]=N*ynew[n]
		#print xnew[n],ynew[n],entropy_L[m][int(N*xnew[n])]*1.0/N	
						

exec "A=np.loadtxt('entropy_2para_8.dat')" 
for i in range (0,len(A)):
	if(entropy_L[int((A[i][1]+N)/2)][int(A[i][0])]==1):
    		entropy_L[int((A[i][1]+N)/2)][int(A[i][0])]=float(A[i][2])
for n in range (0,3*N+1):
	for m in range (0,N+1):
		if(entropy_L[m][n]!=0 and entropy_L[m][n]!=1):
			break
			print m,n,entropy_L[m][n]

for m in range (50,N-50):
	A=[]
	B=[]
	for n in range (0,3*N+1):
		if (entropy_L[m][n]!=1 and entropy_L[m][n]!=0):
			A.append(entropy_L[m][n])
			B.append(n)		
	xnew=[]
	X=np.zeros(shape=len(A))
	Y=np.zeros(shape=len(A))
	for i in range (len(A)):
		X[i]=B[i]*1.0
		Y[i]=A[i]*1.0
 	popt,pcov = optimize.curve_fit(func2,X,Y)	
	flag=0
	for n in range (40,3*N,2):
		if(entropy_L[m][n]==1):
			if(flag==0):
			#	print 2*m-N,n
				delta=entropy_L[m][n-2]-func2(n,popt[0],popt[1],popt[2])
				flag=1
			entropy_L[m][n]=func2(n,popt[0],popt[1],popt[2])+delta
		else:
			flag=0
	
			
for n in range (0,3*N+1):
	for m in range (0,N+1):
		if(entropy_L[m][n]!=0 and entropy_L[m][n]!=1):
			#break
			print 2*m-N,n,entropy_L[m][n]
