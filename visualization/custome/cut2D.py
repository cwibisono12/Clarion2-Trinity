#!/usr/bin/env python3

import matplotlib.pyplot as plt

class cut2D:
	'''
	Class for making 2D cut for a given PID
	Usage:
	1). To Draw a banana gate
	2). To View stored banana gate
	C. Wibisono
	10/23 '23
	'''
	def __init__(self,ax,*,bangate=None):
		'''
		Constructor Description:
		1). ax (required): the axes from a given PID figure
		2). bangate file (optional): the banana gate files, use this argument if we already have
		    banana gate files. 
		    See Clarion2-Trinity/param/2d_all_p.banx as an example of banana gate file
		'''
		self.ax=ax
		self.xcoords=[]
		self.ycoords=[]
		self.index=0
		if bangate !=None:
			banfile=open(bangate)
			lines=banfile.readlines()
			self.lines=[]
			self.idlines=[]
			for line in lines:
				if line.find('#') == -1:
					self.lines.append(line.split())
				else:
					self.idlines.append(line.split())	
			

	def onclick(self,event):
		'''
		Method's Description:
		To draw a banana gate by clicking onto a corresponding PID.
		'''
		if event.inaxes != self.ax.axes:
			return
		self.index=self.index+1
		x,y=event.xdata,event.ydata
		self.xcoords.append(x)
		self.ycoords.append(y)
		self.ax.plot(self.xcoords,self.ycoords,'ro-')
		print(self.index,x,y)
		self.ax.figure.canvas.draw()		

	def cutfile(self,banid):
		'''
		Method's Description:
		To view banana cut files onto a corresponding PID.
		'''
		coordsx=[]
		coordsy=[]
		for i in range(len(self.lines)):
			if int(self.lines[i][0]) == banid:
				coordsx.append(float(self.lines[i][2]))
				coordsy.append(float(self.lines[i][3]))
				self.ax.plot(coordsx,coordsy,'b*-')
							
		self.ax.figure.canvas.draw()

	def connect(self):
		'''
		Method's Description:
		To connect instance onto Matplotlib User Interactive
		'''
		self.ax.figure.canvas.mpl_connect('button_press_event',self.onclick)


	
if __name__ == "__main__":
	import numpy as np
	import sys

	fig,ax=plt.subplots()
	ax.plot(np.arange(0,4096,1),np.arange(0,4096,1))
	banfile=sys.argv[1]	
	banid=int(sys.argv[2])

	p=cut2D(ax,bangate=banfile)
	p.connect()
	p.cutfile(banid)
	plt.show()
