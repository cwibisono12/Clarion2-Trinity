#!/usr/bin/env python3

#Python Fitter Class:

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cvt






class fit:
	def __init__(self,x,y,bkgleft,bkgright,xlow,xup):
		self.x=x
		self.y=y
		self.bkgleft=bkgleft
		self.bkgright=bkgright
		self.xlow=xlow
		self.xup=xup

	def gaussfunc(x,A,mu,sigma,background):
		return A*np.exp(-((x-mu)**2.)/(2*(sigma**2.)))+background

	
	def linear(x,A,B):
		return A*x+B

	def gauss(self):
		backgheight=np.average(y[int(self.bkgleft):int(self.bkgright)])
		mu0=self.x[int(self.xlow)]+(x[int(self.xup)]-x[int(self.xlow)])/2
		meanfitcenter=int(int(self.xlow)+(int(self.xup)-int(self.xlow))/2)
		width=2*(self.x[int(self.xup)]-self.x[meanfitcenter])
		p0=np.array([backgheight,mu0,width])
		popt,pcov=cvt(lambda x, A, mu, sigma: fit.gaussfunc(x,A,mu,sigma,background),self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		return popt
	
	def linfit(self):
		#p0=np.array([2., 0.])
		p0=2.
		#popt,pcov=cvt(fit.linear,self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		popt,pcov=cvt(lambda x, A:fit.linear(x,A,0),self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		return popt




'''
z1=np.arange(0,100,1)
z2=2*z1


#Instantiate the Object:
p=fit(z1,z2,10,20,25,35)

print(p.linfit())
		 
'''
def func():
	print("This is a Class for Fitting")		
