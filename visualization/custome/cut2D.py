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
		self.linelist=[]
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
		lineplot,=self.ax.plot(self.xcoords,self.ycoords,'ro-')
		self.linelist.append(lineplot)
		print(self.index,x,y)
		self.ax.figure.canvas.draw()		

	def onpress(self,event):
		'''
		Method's Description:
		To clear the most recent drawed banana cuts.
		On the axes for a given figure press the button d to delete the most recent 
		banana cut that is just drawn.
		'''
		if event.inaxes != self.ax.axes:
			return
		if event.key != 'd':
			return

		del self.xcoords
		del self.ycoords
		self.xcoords=[]
		self.ycoords=[]
		self.index=0
		for i in range(len(self.linelist)):
			self.linelist.pop(0).remove()

		del self.linelist
		self.linelist=[]
		self.ax.figure.canvas.draw()

	def printbangate(self,event):
		'''
		Method's Description:
		To print coordinates of drawed banana gates onto the screen.
		On the axes for a given figure press the button y to print the most recent banana cuts coordinates.
		'''
		if event.inaxes != self.ax.axes:
			return
		if event.key != 'y':
			return
		print("Printing Drawed Banana Gates onto the screen:\n")
		print("Input Banana ID:\n")
		banid=input()
		print(banid,len(self.xcoords),'#')
		for i in range(len(self.xcoords)):
			print(banid,len(self.xcoords),self.xcoords[i],self.ycoords[i])
			

	def cutfile(self,banid):
		'''
		Method's Description:
		To view banana cut files onto a corresponding PID.
		'''
		self.coordsx=[]
		self.coordsy=[]
		self.lineban=[]
		for i in range(len(self.lines)):
			if int(self.lines[i][0]) == banid:
				self.coordsx.append(float(self.lines[i][2]))
				self.coordsy.append(float(self.lines[i][3]))
				lineplot,=self.ax.plot(self.coordsx,self.coordsy,'b*-')
				self.lineban.append(lineplot)		
		self.ax.figure.canvas.draw()
		
	def clrcutfile(self,event):
		'''
		Method's Description:
		To clear the loaded banana cuts on the given axes on a given figure.
		Press the button k to clear the loaded banana cuts.
		'''
		if event.inaxes !=self.ax.axes:
			return
		if event.key != 'k':
			return

		for i in range(len(self.lineban)):
			self.lineban.pop(0).remove()
		
		del self.lineban
		del self.coordsx
		del self.coordsy
		self.ax.figure.canvas.draw()

	def connect(self):
		'''
		Method's Description:
		To connect instance onto Matplotlib User Interactive
		'''
		self.ax.figure.canvas.mpl_connect('button_press_event',self.onclick)
		self.ax.figure.canvas.mpl_connect('key_press_event',self.onpress)
		self.ax.figure.canvas.mpl_connect('key_press_event',self.printbangate)	
		self.ax.figure.canvas.mpl_connect('key_press_event',self.clrcutfile)	

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
