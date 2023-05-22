#!/usr/bin/env python3

#Python Fitter Class with Interactive Extension:
#User is able to make a Gauss Fit with flat background using mouse click:
#How to Use:
#In order to fit a peak, a plot object has to be instantiated. Example:
#Let p be an object that is instantiated from a given plot characterized by a line object
#1). p=fitgui(line)
#A mouse click event need to be passed from matplotlib base class interactive figure as below:
#2). fig.canvas.mpl_connect('button_press_event',p.onclick)
#A keyboard press event need to be passed as well:
#3). fig.canvas.mpl_connect('key_press_event',p.gauss)
#Interactive Plot is now ready to be used
#How to fit a peak:
#0). Navigate the cursor to the axes where the peak of interest is
#Determining the background region:
#1). Click the left regime of the background
#2). Click the right regime of the background
#Determining the peak region:
#3). Click the left regime of the peak
#4). Click the right regime of the peak
#Fit a peak:
#5). Hit f

#C. Wibisono 05/22/2023

import numpy as np
from matplotlib.lines import Line2D
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
		print(self.xdim,self.ydim)
		self.xlow=0
		self.xup=0
		self.bkgleft=0
		self.bkgright=0
		self.vertical=[]

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

		self.index=self.index+1
		print('index:',self.index,'xdata:',event.xdata,'ydata:',event.ydata,'x:',self.x[int(event.xdata)-min(self.x)],'y(x):',self.y[int(event.xdata)-min(self.x)])
		#print(len(self.vertical))		
		if self.index == 1:
			self.bkgleft = self.x[int(event.xdata)-min(self.x)]
		if self.index == 2:
			self.bkgright = self.x[int(event.xdata)-min(self.x)]
		if self.index == 3:
			self.xlow = self.x[int(event.xdata)-min(self.x)]
		if self.index == 4:
			self.xup = self.x[int(event.xdata)-min(self.x)]
			self.index = 0
			
	def gaussfunc(x,A,mu,sigma,background):
		return A*np.exp(-((x-mu)**2.)/(2*(sigma**2.)))+background
	
	def gauss(self,event):
		if event.inaxes != self.line.axes:
			return
		if event.key != 'f':
			return
		
		for i in range(0,4,1):
			abc2=self.vertical.pop(0)
			abc2.remove()
			#self.vertical[i].remove()
		backgheight=np.average(self.y[int(self.bkgleft)-min(self.x):int(self.bkgright)-min(self.x)])
		#mu0=self.x[int(self.xlow)]+(self.x[int(int(self.xup)-int(self.xlow)/2)])
		meanfitcenter=int(int(self.xlow)+(int(self.xup)-int(self.xlow))/2)
		width=(self.x[int(self.xup)-min(self.x)]-self.x[meanfitcenter-min(self.x)])
		mu0=self.x[meanfitcenter-min(self.x)]
		A0=self.y[meanfitcenter-min(self.x)]-backgheight
		p0=np.array([A0,mu0,width])
		#print(mu0,A0,p0)
		popt,pcov=cvt(lambda x, A, mu, sigma: gui.gaussfunc(x,A,mu,sigma,backgheight),self.x[int(self.xlow)-min(self.x):int(self.xup)-min(self.x)],self.y[int(self.xlow)-min(self.x):int(self.xup)-min(self.x)],p0)
		#return popt,backgheight
		area=popt[0]*np.sqrt(2.0*np.pi)*abs(popt[2])
		print('mean:',popt[1],'sigma:',abs(popt[2]),'FWHM:',abs(popt[2]*gui.FWHMc),'area:',area)
		self.line.axes.plot(self.x[int(self.xlow)-min(self.x):int(self.xup)-min(self.x)],gui.gaussfunc(self.x[int(self.xlow)-min(self.x):int(self.xup)-min(self.x)],*popt,backgheight),'m',linewidth=0.5)
