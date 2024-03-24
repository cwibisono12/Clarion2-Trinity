#!/usr/bin/env python3

import sys
from scipy.signal import find_peaks as fp
import matplotlib.pyplot as plt
import numpy as np

class PeakFinder:
	'''
	Peak Finder class primarily for HPGe calibration.	
	C. Wibisono
	03/24 '24
	'''
	def __init__(self, line):
		self.line = line
		self.x = line.get_xdata()
		self.y = line.get_ydata()
		self.maxheight = max(self.y)
		self.result = None
	
	def connect(self):
		self.line.figure.canvas.mpl_connect('key_press_event',self.find_peak)
	
	def find_peak(self, event):
		if event.inaxes != self.line.axes:
			return
		if event.key != 'h':
			return
		sys.stdout.flush()
		pval,_ = fp(self.y, height = self.maxheight/15., distance = 500)
		if self.result != None:
			self.result.pop(0).remove()

		self.result = self.line.axes.plot(self.x[pval],self.y[pval],"x")
		self.line.figure.canvas.draw()
		n = len(pval)
		for i in range(n):
			print("x:", self.x[pval[i]], "counts:", self.y[pval[i]])


if __name__ == "__main__":
	'''
	Below is an example of how to use Peak Finder Class to determine peaks in Eu152 spectrum
	'''
	from projmod import clarion as p
	from rebin import rebin as r
	from fitgui import fitgui as f

	filemat=sys.argv[1]
	ydim = int(sys.argv[2])
	xdim = int(sys.argv[3])
	xlow = int(sys.argv[4])
	xup = int(sys.argv[5])
	rebinfactor = int(sys.argv[6])

	obj = p(filemat, ydim, xdim, 0, xdim-1, 0, ydim-1)
	objtrp = obj.readparse()
	yobj = obj.project(1,0,1,1,0,0,objtrp)
	yobjrebin = r(xdim, yobj, rebinfactor).rebin()
	fig, ax = plt.subplots()
	x = np.arange(0,xdim,1)
	line, = ax.plot(x[xlow:xup],yobjrebin[xlow:xup],linewidth = 0.85,ls='steps-mid',color='b')
	ax.set_xlim(xlow,xup)
	lineobj = PeakFinder(line)
	lineobj.connect()
	linefit = f(line)
	linefit.connect()
	plt.show()
