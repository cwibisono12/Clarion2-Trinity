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
	#scale=1./(np.power(2.*(np.pi)*(np.power(sigma,2.0)),0.5))
	return np.exp(-(np.power(m,2.0))/(2*np.power(sigma,2.)))

#Normalization of magnetic subtrat:
def normmag(Ji,sigma):
        cm=0.
        for j in range(0,int(4*Ji+2),2):
                mi=0.5*(j-2.*Ji)
                cm=cm+magnetic(mi,sigma)
        return 1/cm

#Bk Coefficient:
def Bk(Ji,k,sigma):
	Bk=0.
	for j in range(0,int(4*Ji+2),2):
		m=0.5*(j-2.*Ji)
		Bk=Bk+normmag(Ji,sigma)*magnetic(m,sigma)*(np.power(-1,Ji-m))*(np.power(2.*Ji+1.,0.5))*CG(Ji,m,Ji,-m,k,0).doit().evalf()
		#print("j:", j,"Bk:", Bk)
	return Bk

#Rk Coefficient:
def Rk(k,l0,l1,Ji,Jf):
	a=(-1.)**(1+(Ji-Jf)+(l1-l0)-k)
	b=np.power(((2.*Ji)+1)*((2.*l0)+1)*((2.*l1)+1),0.5)
	c=CG(l0,1,l1,-1,k,0).doit().evalf()
	d=racah(Ji,Ji,l0,l1,k,Jf).evalf()
	return a*b*c*d

