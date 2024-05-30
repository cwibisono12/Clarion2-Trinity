#!/usr/bin/env python3

#Python Fitter Class with GUI Extension:
#2nd version --Updated on 06/11/2023 --to accomodate compatibility with another plot objects
#C. Wibisono 05/25/2023
#Compatibility:
#This Python Fitter Class with User Interactive extension can be used to fit a peak from any Python object plot.
#How to Use:
#Let line be a curve in a given Python plot.
#1). Instantiate the Object:
#p=fitgui(line)
#2). Connect the object p into Matplotlib User Interactive Binding:
#p.connect()
#3). Ready to fit a peak.

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cvt


class fitgui:
	FWHMc = 2.355

	def __init__(self,line):
		self.index = 0
		self.line = line	
		self.x = line.get_xdata()
		self.y = line.get_ydata()
		self.xdim = len(self.x)
		self.ydim = len(self.y) 
		#print(self.xdim,self.ydim)
		self.xlow=0
		self.xup=0
		self.bkgleft=0
		self.bkgright=0
		self.vertical=[]
		self.gfit = None
		self.curv1 = None
		self.curv2 = None
		self.curv3 = None
		self.tripg1 = None
		self.tripg2 = None
		self.tripg3 = None
		self.tripg4 = None

	def connect(self):
		self.line.figure.canvas.mpl_connect('button_press_event',self.onclick)
		self.line.figure.canvas.mpl_connect('key_press_event',self.gauss)
		self.line.figure.canvas.mpl_connect('key_press_event',self.doubgauss)
		self.line.figure.canvas.mpl_connect('key_press_event',self.tripgauss)

	def print(self):
		print(list(zip(self.x,self.y)))	
	
	def onclick(self,event):
		if event.inaxes != self.line.axes:
			return

		
		vline = self.line.axes.axvline(x=event.xdata)
		self.vertical.append(vline)
		
		if len(self.vertical) > 4:
			for i in range(0,4,1):
				#self.vertical[i].remove()
				#del self.vertical[i]
				abc=self.vertical.pop(0)
				abc.remove()
			#for i in range(0,4,1):
			#	del self.vertical[i]
			self.index=0

		self.line.figure.canvas.draw()		
		self.index=self.index+1
		#print('index:',self.index,'xdata:',event.xdata,'ydata:',event.ydata,'x:',self.x[int(event.xdata)-min(int(self.x))],'y(x):',self.y[int(event.xdata)-min(int(self.x))])
		q=np.where(self.x==round(event.xdata,2))[0][0]   
    #yq=self.y[q]
		print('index:',self.index,'xdata:',event.xdata,'ydata:',event.ydata,'x:',self.x[q],'y(x):',self.y[q])
    #print('index:',self.index,'xdata:',event.xdata,'ydata:',event.ydata,'x:',q,'y(x):',self.y[q])
		#print(len(self.vertical))		
		if self.index == 1:
			#self.bkgleft = self.x[int(event.xdata)-min(self.x)]
			self.bkgleft = round(event.xdata,2)
		if self.index == 2:
			#self.bkgright = self.x[int(event.xdata)-min(self.x)]
			self.bkgright = round(event.xdata,2)
		if self.index == 3:
			#self.xlow = self.x[int(event.xdata)-min(self.x)]
			self.xlow = round(event.xdata,2)
		if self.index == 4:
			#self.xup = self.x[int(event.xdata)-min(self.x)]
			self.xup = round(event.xdata,2)
			self.index = 0

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
	
	def gauss(self,event):
		if event.inaxes != self.line.axes:
			return
		if event.key != 'm':
			return		
		sys.stdout.flush()
		if self.gfit != None:
			abc3=self.gfit.pop(0)
			abc3.remove()
		#backgheight=np.average(self.y[int(self.bkgleft)-min(self.x):int(self.bkgright)-min(self.x)])
		backgheight=np.average(self.y[np.where(self.x==self.bkgleft)[0][0]:np.where(self.x==self.bkgright)[0][0]])
		#mu0=self.x[int(self.xlow)]+(self.x[int(int(self.xup)-int(self.xlow)/2)])
		meanfitcenter=round((self.xlow+(self.xup-self.xlow)/2),2)
		#print('mean:',meanfitcenter)
		z=np.where(self.x==meanfitcenter)[0][0]
    		#width=(self.x[int(self.xup)-min(self.x)]-self.x[meanfitcenter-min(self.x)])
		width=(self.x[np.where(self.x==self.xup)[0][0]]-self.x[z])
		mu0=self.x[z]
		A0=self.y[z]-backgheight
		p0=np.array([A0,mu0,width])
		#print(mu0,A0,p0)
		#popt,pcov=cvt(lambda x, A, mu, sigma: fitgui.gaussfunc(x,A,mu,sigma,backgheight),self.x[int(self.xlow)-min(self.x):int(self.xup)-min(self.x)],self.y[int(self.xlow)-min(self.x):int(self.xup)-min(self.x)],p0)
		popt,pcov=cvt(lambda x, A, mu, sigma: fitgui.gaussfunc(x,A,mu,sigma,backgheight),self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],self.y[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],p0)
		#return popt,backgheight
		area=popt[0]*np.sqrt(2.0*np.pi)*abs(popt[2])
		print('axes:',self.line.axes,'mean:',popt[1],'sigma:',abs(popt[2]),'FWHM:',abs(popt[2]*fitgui.FWHMc),'area:',area)
		self.gfit=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussfunc(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],*popt,backgheight),'m',linewidth=1.5)
		for i in range(0,4,1):
			abc2=self.vertical.pop(0)
			abc2.remove()
		self.line.figure.canvas.draw()
	
	def doubgauss(self,event):	
		if event.inaxes != self.line.axes:
			return
		if event.key != 'd':
			return
		backgheight=np.average(self.y[np.where(self.x==self.bkgleft)[0][0]:np.where(self.x==self.bkgright)[0][0]])	
		print('Enter the 1st width:\n')
		w1=input()
		print('Enter the 2nd width:\n')
		w2=input()
		meanfit1=round((self.xlow+int(int(w1)/2.)),2)
		meanfit2=round((self.xup-int(int(w2)/2.)),2)
		z1=np.where(self.x==meanfit1)[0][0]
		z2=np.where(self.x==meanfit2)[0][0]
		amp1=self.y[z1]
		amp2=self.y[z2]
		sigma1=int(int(w1)/2.)
		sigma2=int(int(w2)/2.)
		p0=np.array([amp1,meanfit1,sigma1,amp2,meanfit2,sigma2])
		popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2: fitgui.gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight),self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],self.y[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],p0)
		#return popt,backgheight	
		area1=popt[0]*np.sqrt(2.0*np.pi)*abs(popt[2])
		area2=popt[3]*np.sqrt(2.0*np.pi)*abs(popt[5])
		areasum=area1+area2
		print('mean1:',popt[1],'sigma1:',abs(popt[2]),'FWHM1:',abs(popt[2]*fitgui.FWHMc),'area1:',area1,'mean2:',popt[4],'sigma2:',abs(popt[5]),'FWHM2:',abs(popt[5]*fitgui.FWHMc),'area2:',area2,'areasum:',areasum)
		
		for i in range(0,4,1):
			abc2=self.vertical.pop(0)
			abc2.remove()
		if self.curv1 != None:
			abc4=self.curv1.pop(0)
			abc4.remove()
		if self.curv2 != None:
			abc5=self.curv2.pop(0)
			abc5.remove()
		if self.curv3 != None:
			abc6=self.curv3.pop(0)
			abc6.remove()
		self.curv1=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussdoub(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],*popt,backgheight),'m',linewidth=1.5)
		self.curv2=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussfunc(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],popt[0],popt[1],popt[2],backgheight),'k',linewidth=1.0)
		self.curv3=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussfunc(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],popt[3],popt[4],popt[5],backgheight),'k',linewidth=1.0)
		self.line.figure.canvas.draw()
	
	def tripgauss(self,event):	
		if event.inaxes != self.line.axes:
			return
		if event.key != 't':
			return
		backgheight=np.average(self.y[np.where(self.x==self.bkgleft)[0][0]:np.where(self.x==self.bkgright)[0][0]]) 	
		print('Enter the 1st width:\n')
		w1=input()
		print('Enter the 2nd width:\n')
		w2=input()
		print('Enter the 3rd width:\n')
		w3=input()
		meanfit1=round((self.xlow+int(int(w1)/2.)),2)
		meanfit2=round((self.xlow+int(int(w1)/2.)+int(int(w2)/2.)),2)
		meanfit3=round((self.xup-int(int(w3)/2.)),2)
		z1=np.where(self.x==meanfit1)[0][0]
		z2=np.where(self.x==meanfit2)[0][0]
		z3=np.where(self.x==meanfit3)[0][0]
		amp1=self.y[z1]
		amp2=self.y[z2]
		amp3=self.y[z3]
		sigma1=int(int(w1)/2.)
		sigma2=int(int(w2)/2.)
		sigma3=int(int(w3)/2.)
		p0=np.array([amp1,meanfit1,sigma1,amp2,meanfit2,sigma2,amp3,meanfit3,sigma3])
		popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3: fitgui.gausstrip(x,A1,mu1,sigma1,A2,mu2,sigma2,A3,mu3,sigma3,backgheight),self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],self.y[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],p0)
		area1=popt[0]*np.sqrt(2.0*np.pi)*abs(popt[2])
		area2=popt[3]*np.sqrt(2.0*np.pi)*abs(popt[5])
		area3=popt[6]*np.sqrt(2.0*np.pi)*abs(popt[8])
		areasum=area1+area2
		#return popt,backgheight	
		print('mean1:',popt[1],'sigma1:',abs(popt[2]),'FWHM1:',abs(popt[2]*fitgui.FWHMc),'area1:',area1)
		print('mean2:',popt[4],'sigma2:',abs(popt[5]),'FWHM2:',abs(popt[5]*fitgui.FWHMc),'area2:',area2)
		print('mean3:',popt[7],'sigma3:',abs(popt[8]),'FWHM2:',abs(popt[8]*fitgui.FWHMc),'area3:',area3)
		
		for i in range(0,4,1):
			abc2=self.vertical.pop(0)
			abc2.remove()
		
		if self.tripg1 != None:
			abc3=self.tripg1.pop(0)
			abc3.remove()
		if self.tripg2 != None:
			abc4=self.tripg2.pop(0)
			abc4.remove()
		if self.tripg3 != None:
			abc5=self.tripg3.pop(0)
			abc5.remove()
		if self.tripg4 != None:
			abc6=self.tripg4.pop(0)
			abc6.remove()
		
		self.tripg1=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gausstrip(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],*popt,backgheight),'m',linewidth=1.5)
		self.tripg2=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussfunc(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],popt[0],popt[1],popt[2],backgheight),'k',linewidth=1.0)
		self.tripg3=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussfunc(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],popt[3],popt[4],popt[5],backgheight),'k',linewidth=1.0)
		self.tripg4=self.line.axes.plot(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],fitgui.gaussfunc(self.x[np.where(self.x==self.xlow)[0][0]:np.where(self.x==self.xup)[0][0]],popt[6],popt[7],popt[8],backgheight),'k',linewidth=1.0)
		self.line.figure.canvas.draw()
	
