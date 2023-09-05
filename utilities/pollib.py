#!/usr/bin/env python3
#Module for Theoretical Polarization Functions
#C. Wibisono 08/31 '23

import numpy as np
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import racah
from numpy.polynomial.legendre import legval 

def Fk(k,Jf,L1,L2,Ji):
	"""
	Fk coefficient based on eq.(4)
	T. Aoki et.al, Atomic Data and Nuclear Data Tables
	23, 349-404 (1979)
	"""
	a=(-1.)**(Jf-Ji-1)
	b=np.power(((2*L1+1)*(2*L2+1)*(2.*Ji+1)),0.5)
	c=CG(L1,1,L2,-1,k,0).doit().evalf()
	d=racah(Ji,Ji,L1,L2,k,Jf).doit().evalf()
	return a*b*c*d

def kappa(k,L1,L2):
	"""
	kappa_nu(l1,l2) coefficient based on
	Table II(b). Fagg and Hanna
	Rev. Mod. Physics. (31) 3 , 1959
	"""
	result = 0.
	if k == 2:
		if L1 == 1 and L2 == 1:
			result=-1./2.
		if L1 == 1 and L2 == 2:
			result=-1./6.
		if L1 == 1 and L2 == 3:
			result=-1./12.	
		if L1 == 2 and L2 == 2:
			result=1./2.
		if L1 == 2 and L2 == 3:
			result=-1./4.
		if L1 == 3 and L2 == 3:
			result=1./3.
	if k == 3:
		if L1 == 1 and L2 == 2:
			result=-1./6.
		if L1 == 1 and L2 == 3:
			result=-1./12.	
		if L1 == 2 and L2 == 2:
			result=0.
		if L1 == 2 and L2 == 3:
			result=1./4.
		if L1 == 3 and L2 == 3:
			result=0.
	if k == 4:
		if L1 == 1 and L2 == 3:
			result=-1./12.	
		if L1 == 2 and L2 == 2:
			result=-1./12.
		if L1 == 2 and L2 == 3:
			result=-1./60.
		if L1 == 3 and L2 == 3:
			result=1./3.
	if k == 5:
		if L1 == 2 and L2 == 3:
			result=-1./20.
		if L1 == 3 and L2 == 3:
			result=0.
	if k == 8:
		if L1 == 3 and L2 == 3:
			result=-1./30.
	return result


def Hk(k,L1,L2,delta,Ji,Jf):
	"""
	Hk(L1,L2) coefficient based on eqs. (1),(2), and (3)
	T. Aoki et.al, Atomic data and Nuclear Data Tables 23, 
	349-404 (1979)
	""" 
	a1=(-1.*kappa(k,L1,L1)*Fk(k,Jf,L1,L1,Ji))
	a2=(2.*delta*kappa(k,L1,L2)*Fk(k,Jf,L1,L2,Ji))
	a3=np.power(delta,2.0)*kappa(k,L2,L2)*Fk(k,Jf,L2,L2,Ji)
	a4=1.+np.power(delta,2.0)
	ak=(a1+a2+a3)/a4
	A1=Fk(k,Jf,L1,L1,Ji)
	A2=2.*delta*Fk(k,Jf,L1,L2,Ji)
	A3=np.power(delta,2.0)*Fk(k,Jf,L2,L2,Ji)
	Ak=(A1+A2+A3)/a4
	Hkcoeff = 2.*ak/Ak
	if Ak == 0 and ak == 0:
		#Hkcoeff = -2*kappa(k,L1,L1)
		Hkcoeff = -1./6.
		#print('Fk:', Fk(k,Jf,L1,L1,Ji))
	print('Hkcoef:\n')
	print('k:', k,'L1:', L1,'L2:', L2,'kappa:', kappa(k,L1,L2),'Hk:', Hkcoeff,'ak:', ak, 'Ak:', Ak)
	return Hkcoeff


def W(phi,theta,a2th,a4th,delta,Ji,Jf,mult):
	"""
	Theoretical polarization distribution based on eqs. (8) 
	T. Aoki et.al, Atomic data and Nuclear Data Tables 23,
	349-404 (1979)
	"""
	#L1=int(abs(Ji-Jf))
	#L2=L1+1
	if abs(round(Ji,1)-round(Jf,1)) > 0:
		L1 = int(abs(Ji-Jf))
		L2 = L1 + 1
	else:
		L1 = 1
		L2 = 2
		
	a0=1.
	a2=a2th*legval(np.cos(np.deg2rad(theta)),[0,0,1,0,0])
	a4=a4th*legval(np.cos(np.deg2rad(theta)),[0,0,0,0,1])
	P22=3.*(1-np.power(np.cos(np.deg2rad(theta)),2.0)) #Associated Leg.Poly for order 2 l=2
	P42=(15./2.)*(7.*np.power(np.cos(np.deg2rad(theta)),2.0)-1.)*(1.-np.power(np.cos(np.deg2rad(theta)),2.0)) #Associated Leg.Polt for order 2 l=4
	A20=Hk(2,L1,L2,delta,Ji,Jf)*P22*a2th
	A40=Hk(4,L1,L2,delta,Ji,Jf)*P42*a4th
	if mult == 1: #Electric L2 type
		phase = 1.
	else:
		phase = -1.
	A=phase*np.cos(np.deg2rad(2.*phi))*0.5*(A20+A40)
	return a0+a2+a4+A


def polpredict(theta,a2th,a4th,delta,Ji,Jf,mult):
	"""
	Theoretical Polarization as a function of scattering angle theta
	P(theta)
	T. Aoki et.al, Atomic data and Nuclear Data Tables 23,
	349-404 (1979)
	"""
	perp=W(0,theta,a2th,a4th,delta,Ji,Jf,mult)
	para=W(90,theta,a2th,a4th,delta,Ji,Jf,mult)
	num=perp-para
	denum=perp+para
	#print('perp:', perp,'para:', para)
	return num/denum


if __name__ == "__main__":
	import sys
	theta=float(sys.argv[1])
	a2th=float(sys.argv[2])
	a4th=float(sys.argv[3])
	delta=float(sys.argv[4])
	Ji=float(sys.argv[5])
	Jf=float(sys.argv[6])
	mult=int(sys.argv[7])
	pol=polpredict(theta,a2th,a4th,delta,Ji,Jf,mult)
	print('theta:', theta,'P(theta):', pol)
