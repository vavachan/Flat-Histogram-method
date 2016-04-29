import numpy as np
import sys
L=int(sys.argv[2])
N=L*L*L
ISING=np.zeros(shape=(L,L,L))
CUBE=np.zeros(shape=N)
N_plus=int(sys.argv[1])
sum=0
sides=np.zeros(shape=3)
a=np.zeros(shape=3)
	
N_plus_t=N_plus
if N_plus>N/2:
	N_plus=N-N_plus
if N_plus==0:
	print (2*N_plus_t-N),(3*N-0),(3*N-0) 
	sys.exit(0)	
def reset():
	for i in range (0,L):
		for j in range (0,L):
			for k in range (0,L):
				ISING[i][j][k]=-1
def energy_n():
	Energy=0
	for i in range (0,L):
		for j in range (0,L):
			for k in range (0,L):
				Ni_plus=i+1
				Ni_minus=i-1
				Nj_plus=j+1
				Nj_minus=j-1
				Nk_plus=k+1
				Nk_minus=k-1
				if (i==0):
				    Ni_minus=L-1
				if (i==L-1):
				    Ni_plus=0
				if (j==0):
				    Nj_minus=L-1
				if (j==L-1):
				    Nj_plus=0
				if (k==0):
				    Nk_minus=L-1
				if (k==L-1):
					Nk_plus=0
#Energy is the -1*sum_{allpaires}(s_i*s_j), hence this function returns the nnumber of unlike bonds.
				Energy=Energy+-1*(ISING[i][j][k]*ISING[i][Nj_plus][k]+ISING[i][j][k]*ISING[i][Nj_minus][k]+ISING[i][j][k]*ISING[Ni_plus][j][k]+ISING[i][j][k]*ISING[Ni_minus][j][k]+ISING[i][j][k]*ISING[i][j][Nk_plus]+ISING[i][j][k]*ISING[i][j][Nk_minus])
	return (3*N+Energy/2)/2

for i in range (0,N):
	if i*i*i>N_plus:
		sum=(i-1)**3
		sides[0]=sides[1]=sides[2]=i-1
		break
while 1:
	a[0]=sides[1]*sides[2]
	a[1]=sides[0]*sides[2]
	a[2]=sides[1]*sides[0]
	ma=max(a)
	vol=sides[0]*sides[1]*sides[2]
	for i in range (0,3):
		if a[i]==ma:
			loca=i
			break
	if vol+a[loca]>N_plus:
		break
	else :
		sides[loca]=sides[loca]+1
vol=sides[0]*sides[1]*sides[2]
rem=N_plus-vol

#print sides
reset()	

for i in range (0,int(sides[0])):
	for j in range (0,int(sides[1])):
		for  k in range (0,int(sides[2])):
			ISING[i][j][k]=1
#print loca,(loca+1)%3,(loca+2)%3
count=0
for i in range (0,int(max(sides[(loca+1)%3],sides[(loca+2)%3]))):
	for j in range (0,int(min(sides[(loca+1)%3],sides[(loca+2)%3]))):
		if(count==rem):
			break
		if loca==0:
			if min(sides[(loca+1)%3],sides[(loca+2)%3])==sides[1]:
				ISING[int(sides[loca])][j][i]=1
			else:
				ISING[int(sides[loca])][i][j]=1
		if loca==1:
			if min(sides[(loca+1)%3],sides[(loca+2)%3])==sides[0]:
				ISING[j][int(sides[loca])][i]=1
			else:
				ISING[i][int(sides[loca])][j]=1
		if loca==2:
			if min(sides[(loca+1)%3],sides[(loca+2)%3])==sides[0]:
				ISING[j][i][int(sides[loca])]=1
			else:
				ISING[i][j][int(sides[loca])]=1
		count=count+1
#			break
			#	print i,j,k
E1=energy_n()
##########################################################################
reset()
count=0
for i in range (0,L):
	for j in range (0,L):
		for k in range (0,L):
			if count==N_plus:
				break
			ISING[i][j][k]=1 
			count=count+1
E2=energy_n()
#########################################################################
R=N_plus%6
S=N_plus/6
#print R,S
sides=np.zeros(shape=2)
a=np.zeros(shape=2)
for i in range (0,L):
	if i*i>S:
		sides[0]=sides[1]=i-1
		break
while 1:
	a[0]=sides[1]#*sides[2]
	a[1]=sides[0]#*sides[2]
	ma=max(a)
	vol=sides[0]*sides[1]
	for i in range (0,2):
		if a[i]==ma:
			loca=i
			break
	if vol+a[loca]>S:
		break
	else :
		sides[loca]=sides[loca]+1
rem=S-vol
#print sides
reset()
for i in range (0,int(sides[0])):
	for j in range (0,int(sides[1])):
		for  k in range (0,L):
			ISING[i][j][k]=1
#print ISING
#print '******'
if(sides[0]>sides[1]):
	for i in range (0,int(rem)):
		for  k in range (0,L):
			ISING[i][sides[1]][k]=1
	for k in range  (0,R):
		ISING[rem][sides[1]][k]=1	
else :
	for j in range (0,int(rem)):
		for  k in range (0,L):
			ISING[sides[0]][j][k]=1
	for k in range  (0,R):
		ISING[sides[0]][rem][k]=1	
E3=energy_n()
##########################################################################
#print ISING
reset()
sq=L*L
S=N_plus_t/sq
R=N_plus_t%sq
#print S,R
for i in range (0,int(S)):
	for j in range (0,L):
		for k in range (0,L):
			ISING[i][j][k]=1	
sides=np.zeros(shape=2)
a=np.zeros(shape=2)
for i in range (0,L):
	if i*i>R:
		sides[0]=sides[1]=i-1
		break
while 1:
	a[0]=sides[1]#*sides[2]
	a[1]=sides[0]#*sides[2]
	ma=max(a)
	vol=sides[0]*sides[1]
	for i in range (0,2):
		if a[i]==ma:
			loca=i
			break
	if vol+a[loca]>R:
		break
	else :
		sides[loca]=sides[loca]+1
rem=R-vol
for j in range (0,int(sides[0])):
	for k in range (0,int(sides[1])):
		ISING[S][j][k]=1
if(sides[0]>sides[1]):
	for i in range (0,int(rem)):
		ISING[S][i][sides[1]]=1
else:
	for i in range (0,int(rem)):
		ISING[S][sides[0]][i]=1
E4=energy_n()
print(float(2*N_plus_t-N),(3*N-min(E1,E2,E3,E4)),float(3*N-6*N_plus))
#3*N-E1,3*N-E2,3*N-E3,3*N-E4#(3*N-min(E1,E2,E3,E4)),(3*N-min(E1,E2,E3)),float(3*N-6*N_plus)
# the magnetization, number of like bonds. 
