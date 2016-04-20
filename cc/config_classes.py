import numpy as np
import math
L1=6
N1=L1**3
exec "A=np.loadtxt('cc_%d.dat',skiprows=0)" %L1
entropy_L1=np.zeros(shape=(N1+1,3*N1+1))
entropy_L1[0][3*N1]=1
entropy_L1[1][3*N1-6]=1
entropy_L1[N1][3*N1]=1
entropy_L1[N1-1][3*N1-6]=1
for m in range (0,N1+1):
        for n in range (int(A[m][2]),int(A[m][1])+1,2):
                entropy_L1[m][n]=1
for m in range (0,N1+1):
        for n in range (0,3*N1+1):
		if (entropy_L1[m][n]!=0):
			print n,2*m-N1,entropy_L1[m][n]






