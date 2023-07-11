#!/usr/bin/env python3
#Functions to Calculate Bk and Rk coeff
#Useful for Angular Distribution
#C. Wibisono 07/11 '23

import numpy as np
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import racah
from sympy import *
import sys
import os



def magnetic(m,sigma):
	scale=1./(np.power(2.*(np.pi)*(np.power(sigma,2.0)),0.5))
	return scale*np.exp(-(np.power(m,2.0))/(2*np.power(sigma,2.)))

#Bk Coefficient:
def Bk(Ji,k,sigma):
	Bk=0.
	for j in range(-Ji,Ji+1,1):
		Bk=Bk+magnetic(j,sigma)*(np.power(-1,Ji-j))*(np.power(2.*Ji+1.,0.5))*CG(Ji,j,Ji,-j,k,0).doit().evalf()
		#print("j:", j,"Bk:", Bk)
	return Bk

#Rk Coefficient:
def Rk(k,l0,l1,Ji,Jf):
	a=(-1.)**(1+(Ji-Jf)+(l1-l0)-k)
	b=np.power(((2.*Ji)+1)*((2.*l0)+1)*((2.*l1)+1),0.5)
	c=CG(l0,1,l1,-1,k,0).doit().evalf()
	d=racah(Ji,Ji,l0,l1,k,Jf).evalf()
	return a*b*c*d

