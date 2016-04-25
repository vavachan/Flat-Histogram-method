import numpy as np
import sys
L=int(sys.argv[1])
EN=L*L*L
exec "A=np.loadtxt('entropy_2para_%d.dat',skiprows=0)" %L
def mag(t):
	P_fnc=0
	num_1=0
	num_2=0
	for i in range (0,len(A)):
		P_fnc=P_fnc+A[i][3]*np.exp(-1/(t)*(6*EN-2*A[i][0]))
		num_1=num_1+abs(A[i][1])*A[i][3]*np.exp(-1/(t)*(6*EN-2*A[i][0]))
		num_2=num_2+(A[i][1])**2*A[i][3]*np.exp(-1/(t)*(6*EN-2*A[i][0]))
#		print A[i][0],' ',(A[i][1]+EN)/2,' ',num_1#A[i][3]*np.exp(-1/(t)*(6*EN-2*A[i][0])),P_fnc
	return ((num_2/P_fnc)-(num_1/P_fnc)**2)/t,num_1/P_fnc,P_fnc
exec "infile=open('mag1_%d.dat','w')" %L
for i in range (100,700,10):
	x=i/100.0
	print x
	infile.write(str(x)+'\t'+str(mag(x)[0]/EN)+'\t'+str(mag(x)[1]/EN)+'\t'+str(mag(x)[2]/EN)+'\n')
