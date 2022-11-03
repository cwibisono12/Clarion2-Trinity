#!/usr/bin/env python3

#Python Fitter Class:
#C. Wibisono 11/01/2022

import numpy as np
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
	
	def gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,background):
		gauss1=A1*np.exp(-((x-mu1)**2.)/(2*(sigma1**2.)))	 
		gauss2=A2*np.exp(-((x-mu2)**2.)/(2*(sigma2**2.)))
		return gauss1+gauss2+background
	
	def gausstrip(x,A1,mu1,sigma1,A2,mu2,sigma2,A3,mu3,sigma3,background):	
		gauss1=A1*np.exp(-((x-mu1)**2.)/(2*(sigma1**2.)))	 
		gauss2=A2*np.exp(-((x-mu2)**2.)/(2*(sigma2**2.)))
		gauss3=A3*np.exp(-((x-mu3)**2.)/(2*(sigma3**2.)))
		return gauss1+gauss2+gauss3+background


	def linear(x,A,B):
		return A*x+B

	def gauss(self):
		backgheight=np.average(self.y[int(self.bkgleft):int(self.bkgright)])
		#mu0=self.x[int(self.xlow)]+(self.x[int(int(self.xup)-int(self.xlow)/2)])
		meanfitcenter=int(int(self.xlow)+(int(self.xup)-int(self.xlow))/2)
		width=(self.x[int(self.xup)]-self.x[meanfitcenter])
		mu0=self.x[meanfitcenter]
		A0=self.y[meanfitcenter]-backgheight
		p0=np.array([A0,mu0,width])
		popt,pcov=cvt(lambda x, A, mu, sigma: fit.gaussfunc(x,A,mu,sigma,backgheight),self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		return popt,backgheight
	
	def linfit(self):
		#p0=np.array([2., 0.])
		p0=2.
		#popt,pcov=cvt(fit.linear,self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		popt,pcov=cvt(lambda x, A:fit.linear(x,A,0),self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		return popt

	def doubgauss(self,amp1,meanfit1,sigma1,amp2,meanfit2,sigma2):	
		backgheight=np.average(self.y[int(self.bkgleft):int(self.bkgright)])	
		p0=np.array([amp1,meanfit1,sigma1,amp2,meanfit2,sigma2])
		popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2: fit.gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight),self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		return popt,backgheight	
	
	def tripgauss(self,amp1,meanfit1,sigma1,amp2,meanfit2,sigma2,amp3,meanfit3,sigma3):	
		backgheight=np.average(self.y[int(self.bkgleft):int(self.bkgright)])	
		p0=np.array([amp1,meanfit1,sigma1,amp2,meanfit2,sigma2,amp3,meanfit3,sigma3])
		popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3: fit.gausstrip(x,A1,mu1,sigma1,A2,mu2,sigma2,A3,mu3,sigma3,backgheight),self.x[int(self.xlow):int(self.xup)],self.y[int(self.xlow):int(self.xup)],p0)
		return popt,backgheight	
	



