import numpy as np
import math
from scipy import optimize
import scipy
L2=8
L1=6
#L1 is the linear size of the ising model for which we have the entropy function.L2 is the size of the ising model for which we would like to get a
#initial function for the entropy.
N1=L1**3
N2=L2**3
def func(x,a,b,c):
        return a*x**2+b*x+c
def func1(x,a,b):
        return a*x+b
def func2(x,a):
	return a*x
#N1 and N2 are the no of sites
exec "A=np.loadtxt('%d.dat',skiprows=0)" %L2
#loads the maximum and minimum n={no of unlike bonds} for a corresponding m={no of upspins#}. 
#print A
entropy_L2=np.zeros(shape=(N2+1,3*N2+1))
entropy_L1=np.zeros(shape=(N1+3,3*N1+3))
extended_L1=np.zeros(shape=(N1+3,3*N1+3))
#Extended L1 is an intermediate entropy function which we use for fitting.
for m in range (0,N2+1):
	for n in range (int(A[m][3]),int(A[m][1]),2):
		entropy_L2[m][n]=1
#now entropy_L2 has non-zero value if the configuration class is allowed.
entropy_L2[0][3*N2]=1
entropy_L2[1][3*N2-6]=1
entropy_L2[N2][3*N2]=1
entropy_L2[N2-1][3*N2-6]=1
for m in range (0,N2+1):
	for n in range (0,3*N2+1):
		if(entropy_L2[m][n]):
			break
#			print m,n,1
#m is the number of upspins , and n is the number of like bonds.
for m in range (0,N2+1):
	for n in range (0,3*N2+1):
		if(entropy_L2[m][n]):
#			break
		#	print float(m),float(n),1#
#			print float(int(float(m)/N2*N1))/N1,float(int(float(n)/N2*N1))/N1,1
#			print float(int(float(m)/N2*N1)+1)/N1,float(int(float(n)/N2*N1)+1)/N1,1
			if(int(float(n)/N2*N1)%2==0):
				N=int(float(n)/N2*N1)
			else:
				N=int(float(n)/N2*N1)+1
			extended_L1[int(float(m)/N2*N1)+1][N+2]=1
			extended_L1[int(float(m)/N2*N1)][N]=1
			extended_L1[int(float(m)/N2*N1)][N+2]=1
			extended_L1[int(float(m)/N2*N1)+1][N]=1
#extended_L1 is an entropy function which has non-zero values such that the normalized entropy funciton for L2 always had 4 points of extended_L1 surrounding it.
exec "A=np.loadtxt('entropy_mp_9.dat',skiprows=0)" 
#A has now the entropy function corresponding to #L1
#print A
for i in range (0,len(A)):
#	print (1.0*A[i][1]),(1.0*A[i][0]),float(A[i][2])/N1
	entropy_L1[int((A[i][1]+N1)/2)][int(A[i][0])]=float(A[i][2])/N1
#we write them to entropy_L2 
for n in range (0,3*N1+1):
	for m in range (0,N1/2+1):
		mid=(entropy_L1[m][n]+entropy_L1[N1-m][n])/2
		entropy_L1[m][n]=entropy_L1[N1-m][n]=mid
#making it symetric around m=N/2
X1=np.zeros(shape=3*N1+2)
Y1=np.zeros(shape=3*N1+2)
for m in range (0,N1+2):
	X=[]
	E=[]
	for n in range (0,3*N1+2):
		if(entropy_L1[m][n]!=0):
			E.append(entropy_L1[m][n])
			X.append(float(n)/N1)	
	E_1=np.zeros(shape=len(X))
	X_1=np.zeros(shape=len(X))
	if(len(X) >= 3):
		for i in range (0,len(X)):
			E_1[i]=E[i]
			X_1[i]=X[i]	
			#print E_1[i],X_1[i]
		popt,pcov = optimize.curve_fit(func,X_1,E_1)
	#	for i in range (0,len(X)):
	#		print X_1[i],E_1[i],func(X_1[i],popt[0],popt[1],popt[2])
	########for i in range (0,3*N1+2,1):
	########	if(extended_L1[m][i]!=0):
	########		y=i*1.0/N1
	########	#	X1[i]=y
	########	#	Y1[i]=func(y,popt[0],popt[1],popt[2])
	########		#print float(n),2*float(m)-N1,extended_L1[m][i]
	########		extended_L1[m][i]=func(y,popt[0],popt[1],popt[2])
	if (len(X)==2) :
		for i in range (0,len(X)):
			E_1[i]=E[i]
			X_1[i]=X[i]	
			#print E_1[i],X_1[i]
		popt,pcov = optimize.curve_fit(func1,X_1,E_1)
	if(len(X)==1):
		for i in range (0,len(X)):
			E_1[i]=E[i]
			X_1[i]=X[i]	
			#print E_1[i],X_1[i]
		popt,pcov = optimize.curve_fit(func2,X_1,E_1)


	for i in range (0,3*N1+2,1):
		if(extended_L1[m][i]!=0):
			y=i*1.0/N1
		#	X1[i]=y
		#	Y1[i]=func(y,popt[0],popt[1],popt[2])
			#print float(n),2*float(m)-N1,extended_L1[m][i]
			if (len(X)>=3):
				extended_L1[m][i]=func(y,popt[0],popt[1],popt[2])	
			if (len(X)==2) :
				extended_L1[m][i]=func1(y,popt[0],popt[1])
			if (len(X)==1) :
				extended_L1[m][i]=func2(y,popt[0])
	
			
for m in range (0,N1+2):
	for n in range (0,3*N1+2):
		if(extended_L1[m][n]):
			break
			print float(n),2*float(m)-N1,extended_L1[m][n]
m=1
n=3*N2-6
x=.25*(extended_L1[int(float(m)/N2*N1)][int(float(n)/N2*N1)]+extended_L1[int(float(m)/N2*N1)+1][int(float(n)/N2*N1)]+extended_L1[int(float(m)/N2*N1)][int(float(n)/N2*N1)+1]+extended_L1[int(float(m)/N2*N1)+1][int(float(n)/N2*N1)+1])	
alpha=math.log(N2)/N2

for n in range (0,3*N2+1):
	for m in range (0,N2+1):
		if(entropy_L2[m][n]):
		#	break
			entropy_L2[m][n]=.25*(extended_L1[int(float(m)/N2*N1)][int(float(n)/N2*N1)]+extended_L1[int(float(m)/N2*N1)+1][int(float(n)/N2*N1)]+extended_L1[int(float(m)/N2*N1)][int(float(n)/N2*N1)+1]+extended_L1[int(float(m)/N2*N1)+1][int(float(n)/N2*N1)+1])	
			if(m==0 or m==N2):
				entropy_L2[m][n]=0
			print float(n),2*float(m)-N2,entropy_L2[m][n]#*1.0*alpha/x
			# number of like bonds , magnetization

########for n in range (0,3*N1+1):
########	if(entropy_L1[m][n]):
########		break
########		print float(m)/N1,float(n)/N1,1

