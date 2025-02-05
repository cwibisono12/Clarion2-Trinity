#!/usr/bin/env python3

import numpy as np

'''
Geometrical Attenuation Factor Module
Based on Qk Subroutine Fortran 77 AD Program
Rewritten in Python
02/05 '25
'''

def tau(E):
	EMeV=E/1000.
	ELN=np.log(EMeV)
	TLN=-1.1907-0.5372*ELN-0.043*(ELN**2.)+0.0218*(ELN**3.)+0.0765*(ELN**4.)+0.0095*(ELN**5.)
	return np.exp(TLN)

def terma(R,D,T,E):
	alpha = np.arctan(R/(D+T))
	gamma = np.arctan(R/D)
	deltax1 = alpha/1000.
	t1 = 0
	t2 = 0
	t3 = 0
	for i in range(1001):
		beta = i*deltax1
		if i == 0 or i == 1000:
			a = 1.0
		else:
			if i % 2 == 0:
				a = 2.0
			else:
				a = 4.0

		fac1 = -1.*tau(E)*T/(np.cos(beta))
		temp1 = 0.5*a*(3*((np.cos(beta))**2.)-1.)*(1.-np.exp(fac1))*(np.sin(beta))*deltax1
		temp2 = 0.125*a*(35*((np.cos(beta))**4.)-30.*((np.cos(beta))**2.)+3.)*(1.-np.exp(fac1))*np.sin(beta)*deltax1
		temp3 = a*(1.-np.exp(fac1))*(np.sin(beta))*deltax1

		t1 = t1 +temp1
		t2 = t2 +temp2
		t3 = t3 + temp3

	t1 = t1/3.
	t2 = t2/3.
	t3 = t3/3.

	return [t1, t2, t3]


def termb(R,D,T,E):
	alpha = np.arctan(R/(D+T))
	gamma = np.arctan(R/D)
	deltax2 =(gamma- alpha)/1000.
	t4 = 0
	t5 = 0
	t6 = 0
	for i in range(1001):
		beta = alpha + i*deltax2
		if i == 0 or i == 1000:
			a = 1.0
		else:
			if i % 2 == 0:
				a = 2.0
			else:
				a = 4.0

		fac2 = -1.*tau(E)*((R/np.sin(beta))-(D/np.cos(beta)))
		temp4 = 0.5*a*(3.*((np.cos(beta))**2.)-1.)*(1.-np.exp(fac2))*(np.sin(beta))*deltax2
		temp5 = 0.125*a*(35.*((np.cos(beta))**4.)-30.*((np.cos(beta))**2.)+3.)*(1.-np.exp(fac2))*np.sin(beta)*deltax2
		temp6 = a*(1.-np.exp(fac2))*(np.sin(beta))*deltax2

		t4 = t4 + temp4
		t5 = t5 + temp5
		t6 = t6 + temp6

	t4 = t4/3.
	t5 = t5/3.
	t6 = t6/3.

	return [t4, t5, t6]



if __name__ == "__main__":
	import sys
	import QK
	R = float(sys.argv[1])
	D = float(sys.argv[2])
	T = float(sys.argv[3])
	E = int(sys.argv[4])

	ans1 = terma(R,D,T,E)
	ans2 = termb(R,D,T,E)
	Q2 = (ans1[0]+ans2[0])/(ans1[2]+ans2[2]) #Qk2
	Q4 = (ans1[1]+ans2[1])/(ans1[2]+ans2[2]) #Qk4
	Qkc=QK.qk(E,R,D,T)
	print("Reversed_Engineer:")
	print('Q2:',Q2,'Q4:',Q4)
	print("Fortran:")
	print('Q2:',Qkc[0],'Q4:',Qkc[1])
	
