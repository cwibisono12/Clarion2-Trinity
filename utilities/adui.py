#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

class adui:
	"""
	Usage:
	Class for determining a local minima for a given plot given the lower and upper limit of the domain for a given plot
	C. Wibisono
	08/23 '23
	How to Use:
	1). Instantiate an object p=adui(pythonplot)
	2). Connect to the UI p.connect()
	3). Determine lower and upper region by clicking the mouse button
	4). Hit d to determine local minima
	"""

	def __init__(self,line):
		self.line=line
		self.index=0
		self.vertical=[]
		self.xline=line.get_xdata(orig=True).round(1)
		self.yline=line.get_ydata(orig=True)
		self.lowlim=0
		self.uplim=0
		self.indlowlim=0
		self.induplim=0

	def connect(self):
		self.line.figure.canvas.mpl_connect('button_press_event',self.onclick)
		self.line.figure.canvas.mpl_connect('key_press_event',self.localmin)

	def onclick(self,event):
		if event.inaxes != self.line.axes:
			return
		vline = self.line.axes.axvline(x=event.xdata)
		self.vertical.append(vline)
		if len(self.vertical) > 2:
			for i in range(0,2,1):
				abc=self.vertical.pop(0)
				abc.remove()
			self.index = 0
		self.line.figure.canvas.draw()
		self.index = self.index + 1
		#q = np.where(self.xline == round(event.xdata,1))[0][0]
		q = np.argmin(np.abs(self.xline - round(event.xdata,2)))
		print('index:',self.index,'xdata:',event.xdata,'ydata:',event.ydata,'x:',self.xline[q],'y:',self.yline[q])
		#print('index:',self.index,'xdata:',event.xdata,'ydata:',event.ydata)
		if self.index == 1:
			self.lowlim = round(event.xdata,1)
			self.indlowlim = q
			print('index:', self.indlowlim,'val:', self.lowlim)
			#print(self.xline[100],self.yline[100])
		if self.index == 2:
			self.uplim = round(event.xdata,1)
			self.induplim = q
			print('index:', self.induplim,'val:', self.uplim)
			self.index = 0
		
	def localmin(self,event):
		if event.key != 'd':
			return

		interval=self.induplim-self.indlowlim+1
		#Create new array from the original array but only select portion of it
		newx=np.zeros(interval)
		newy=np.zeros(interval)
		for j in range(0,interval,1):
			newx[j]=self.xline[self.indlowlim+j]
			newy[j]=self.yline[self.indlowlim+j]
		#Finding the local minima index:
		#indlocminima=np.where(newy == np.min(newy))[0][0]
		indlocminima=np.argmin(np.abs(newy - np.min(newy)))
		#The x value asscociates with the local minima:
		xmin=newx[indlocminima]
		print("Local Minima Associated with the Interval\n")
		print('xmin:', xmin,'y(x):', newy[indlocminima])
