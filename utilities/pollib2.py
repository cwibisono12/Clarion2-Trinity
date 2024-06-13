#!/usr/bin/env python3
#Module for Theoretical Polarization
#C. Wibisono 10/02 '23
#This Module does not require AD's coeffs.
#Formulations derived from the following references:
#1).Rev. Mod. Physics (31) 3 1959
#2).Phys. Rev. (88) 4 1952
#3).Satchler, 1955 Polarization in Nuclear Reactions
#4).Definition of Predicted Polarization follows from: Atomic Data and Nuclear Tables (23)
#349 - 404 (1979)

import numpy as np
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import racah
from numpy.polynomial.legendre import legval 
import adlib as ad



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
		if L1 == 3 and L2 == 4:
			result = -1./3.
		if L1 == 4 and L2 == 4:
			result = 5./17.

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
		if L1 == 3 and L2 == 4:
			result = 0.133333
		if L1 == 4 and L2 == 4:
			result = 0
	if k == 4:
		if L1 == 1 and L2 == 3:
			result=-1./12.	
		if L1 == 2 and L2 == 2:
			result=-1./12.
		if L1 == 2 and L2 == 3:
			result=-1./60.
		if L1 == 3 and L2 == 3:
			result=1./3.
		if L1 == 3 and L2 == 4:
			result = -0.022222222
		if L1 == 4 and L2 == 4:
			result = 0.1111111111
	if k == 5:
		if L1 == 2 and L2 == 3:
			result=-1./20.
		if L1 == 3 and L2 == 3:
			result=0.
	if k == 6:
		if L1 == 3 and L2 == 3:
			result=-1./30.
	return result


def W(phi,theta,delta,Ji,Jf,sigma,mult):
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
		
	if mult == 1: #Electric L2 type
		phase = 1.
	else:
		phase = -1.
	
	a0=1.+(delta**2.)
	#a0=1. //revised result from the algebra
	B20=ad.Bk(Ji,2,sigma)
	R20=ad.Rk(2,L1,L1,Ji,Jf)
	R22=2.*ad.Rk(2,L1,L2,Ji,Jf)*delta
	R24=ad.Rk(2,L2,L2,Ji,Jf)*delta**2.
	B40=ad.Bk(Ji,4,sigma)
	R40=ad.Rk(4,L1,L1,Ji,Jf)
	R42=2.*ad.Rk(4,L1,L2,Ji,Jf)*delta
	R44=ad.Rk(4,L2,L2,Ji,Jf)*delta**2.
	#Unpolarized a2 and a4 coefficients:
	a20=B20*(R20+R22+R24)
	a40=B40*(R40+R42+R44)
	#Polarized a2 and a4 coefficients:
	a20pol=phase*B20*(-kappa(2,L1,L1)*R20+kappa(2,L1,L2)*R22+kappa(2,L2,L2)*R24)
	a40pol=phase*B40*(-kappa(4,L1,L1)*R40+kappa(4,L1,L2)*R42+kappa(4,L2,L2)*R44)
	#Unpolarized Term:
	a2=a20*legval(np.cos(np.deg2rad(theta)),[0,0,1,0,0])
	a4=a40*legval(np.cos(np.deg2rad(theta)),[0,0,0,0,1])
	P22=3.*(1-np.power(np.cos(np.deg2rad(theta)),2.0)) #Associated Leg.Poly for order 2 l=2
	P42=(15./2.)*(7.*np.power(np.cos(np.deg2rad(theta)),2.0)-1.)*(1.-np.power(np.cos(np.deg2rad(theta)),2.0)) #Associated Leg.Polt for order 2 l=4
	#Polarized Term:
	Apol=np.cos(np.deg2rad(2.*phi))
	A2=a20pol*P22*Apol
	A4=a40pol*P42*Apol
	Wpol=(a0+a2+a4)+(A2+A4)
	#Wunpol=(a0+a2+a4)
	return Wpol


def polpredict(theta,delta,Ji,Jf,sigma,mult):
	"""
	Theoretical Polarization as a function of scattering angle theta
	P(theta)
	T. Aoki et.al, Atomic data and Nuclear Data Tables 23,
	349-404 (1979)
	"""
	perp=W(0,theta,delta,Ji,Jf,sigma,mult)
	para=W(90,theta,delta,Ji,Jf,sigma,mult)
	num=perp-para
	denum=perp+para
	#print('perp:', perp,'para:', para)
	return num/denum


if __name__ == "__main__":
	import sys
	theta=float(sys.argv[1])
	delta=float(sys.argv[2])
	Ji=float(sys.argv[3])
	Jf=float(sys.argv[4])
	mult=int(sys.argv[5])
	sigma=float(sys.argv[6])
	pol=polpredict(theta,delta,Ji,Jf,sigma,mult)
	print('theta:', theta,'P(theta):', pol)
