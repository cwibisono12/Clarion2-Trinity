#!/usr/bin/env python3

import projmod as p
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tck
import fitgui3 as f
import rebinv2 as r

class gggateplotter(p.clarion):
	'''
	Usage:
	Custome Class to visualize a projection from gg matrix given the gamma matrix and gate file
	C. Wibisono
	10/16/'23
	How to Use:
	1. Object Instantiation:
		Let p be an object:
		p=gggateplotter(ggmat,5000,5000,0,4999,0,4999,gammagatefile)
	2. Displaying the Projected Spectra from gammagatefile:
		p.plot(0,4999,1)
	Ready to fit the peak
	'''	


	plt.rcParams['font.family']='serif'
	plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']	
	plt.rcParams['figure.dpi'] = 200

	def __init__(self,infile,ydim,xdim,xlow,xup,ylow,yup,gatefile):
		'''	
		Constructor Description:
		1). infile ---> ggmatrix file
		2). ydim ----> ydimension of ggmatrix
		3). xdim ---> xdimension of ggmatrix
		4). xlow---> lower x 
		5). xup----> upper x
		6). ylow--->lower y
		7). yup----> upper y
		8). gatefile ---> gamma gate file consisting the region for gating and the region for background
		'''
		
		p.clarion.__init__(self,infile,ydim,xdim,xlow,xup,ylow,yup)
		self.ytrp=self.readparse()
		self.gatefile=open(gatefile)
		gateline=self.gatefile.readlines()
		self.gate=[]
		for line in gateline:
			if line.find('#') == -1:
				self.gate.append(list(line.split()))
			else:
				self.nuclei=list(line.split())			


		self.gatenum=len(self.gate)
		self.x=np.arange(0,xdim,1)
		self.gatefile.close()

	def plot(self,xlow,xup,rebinfactor):
		'''
		Usage:
		Visualizing the projected gamma spectra from the given gatefile and a gamma matrix
		Method's Argument(s):
		xlow: the lower limit of x axis to display the spectra
		xup: the upper limit of y axis to display the spectra
		rebinfactor: the rebinfactor per keV
		'''
		project=[]
		self.projectrebin=[]
		
		for i in range(self.gatenum):
			project.append(self.project(int(self.gate[i][0]),int(self.gate[i][1]),int(self.gate[i][2]),int(self.gate[i][3]),int(self.gate[i][4]),int(self.gate[i][5]),self.ytrp))
			self.projectrebin.append(r.rebin(self.xdim,project[i],rebinfactor).rebin())
		
		self.x = r.rebin(self.xdim,project[0],rebinfactor).rebinx()			

		self.fig,self.ax=plt.subplots(self.gatenum,1)
		self.fig.suptitle(str(self.nuclei[1])+' '+'Gamma Spectra'+'\n'+'Clarion2-Trinity')
		for j in range(self.gatenum):
			if self.gatenum > 1:
				self.ax[j].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				self.ax[j].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				self.ax[j].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.ax[j].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.x = r.rebin(self.xdim, project[0],rebinfactor).rebinx()
			else:
				self.ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				self.ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				self.ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				
		
		self.rebinfactor=rebinfactor	
		self.lineplot=[]
		colorlist=['r','g','b','c','k','m','r','g','b','c','k']
		#Transform xlow and xup in terms of xrebin, and yrebin coordinates space, i.e find the index in xrebin and yrebin given xlow,xup:
		xlowtrf = np.argmin(np.abs(self.x - xlow))
		xuptrf = np.argmin(np.abs(self.x - xup))
		for k in range(0,self.gatenum,1):
			if self.gatenum > 1 :
				lineplotind,=self.ax[k].plot(self.x[xlowtrf:xuptrf],self.projectrebin[k][xlowtrf:xuptrf],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
				self.lineplot.append(lineplotind)
				#lineplot.append(ax[k].plot(self.x[xlow:xup],projectrebinr[k][xlow:xup],linewidth=0.85,ls='steps-mid',color='r',label=str(self.gate[k][6])))
				self.ax[k].legend()
				self.ax[k].set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				self.ax[k].set_ylim(bottom=0)
				self.ax[k].set_xlim(xlow,xup)
			else:
				lineplotind,=self.ax.plot(self.x[xlowtrf:xuptrf],self.projectrebin[k][xlowtrf:xuptrf],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
				self.lineplot.append(lineplotind)
				#lineplot.append(ax[k].plot(self.x[xlow:xup],projectrebinr[k][xlow:xup],linewidth=0.85,ls='steps-mid',color='r',label=str(self.gate[k][6])))
				self.ax.legend()
				self.ax.set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				self.ax.set_ylim(bottom=0)
				self.ax.set_xlim(xlow,xup)
			
		if self.gatenum > 1:
			self.ax[self.gatenum-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
		else:
			self.ax.set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
		
		self.fitline=[]
		for m in range(0,self.gatenum,1):
			self.fitline.append(f.fitgui(self.lineplot[m]))
			self.fitline[m].connect()

		self.fig.canvas.mpl_connect('key_press_event',self.update)
		self.fig.canvas.mpl_connect('key_press_event',self.modify)
		self.fig.canvas.mpl_connect('key_press_event',self.savefigure)
		plt.show()

	def update(self,event):
		'''
		Usage:
		To update plot region based on user's input xlow and xup.
		Press e to expand the plot's region. Enter as follow:
		xlow xup then hit enter/return.
		Example:
		1000 2000
		The new plot's region will be drawn based on the limit provided above.
		'''
		if event.key != 'e':
			return
			
		for i in range(len(self.lineplot)):
			self.lineplot.pop(0).remove()

		del self.lineplot
		print("Input the new xlow and xup\n")
		xlow, xup=map(int,input().split())
	
		#Transform xlow and xup in terms of xrebin and yrebin indices:
		xlowtrf = np.argmin(np.abs(self.x - xlow))
		xuptrf = np.argmin(np.abs(self.x - xup))

		self.lineplot=[]
		colorlist=['r','g','b','c','k','m','r','g','b','c','k']
		for j in range(self.gatenum):
			if self.gatenum > 1:
				self.ax[j].clear()
				self.ax[j].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				self.ax[j].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				self.ax[j].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.ax[j].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
			else:
				self.ax.clear()
				self.ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				self.ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				self.ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

		for k in range(0,self.gatenum,1):
			if self.gatenum > 1:
				lineplotind,=self.ax[k].plot(self.x[xlowtrf:xuptrf],self.projectrebin[k][xlowtrf:xuptrf],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
				self.lineplot.append(lineplotind)
				self.ax[k].legend()
				self.ax[k].set_ylabel('counts/'+str(self.rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				self.ax[k].set_ylim(bottom=0)
				self.ax[k].set_xlim(xlow,xup)
			else:
				lineplotind,=self.ax.plot(self.x[xlowtrf:xuptrf],self.projectrebin[k][xlowtrf:xuptrf],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
				self.lineplot.append(lineplotind)
				self.ax.legend()
				self.ax.set_ylabel('counts/'+str(self.rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				self.ax.set_ylim(bottom=0)
				self.ax.set_xlim(xlow,xup)
	
		if self.gatenum > 1:
			self.ax[self.gatenum-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
		else:
			self.ax.set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
			

		del self.fitline
		self.fitline=[]
		for m in range(0,self.gatenum,1):
			self.fitline.append(f.fitgui(self.lineplot[m]))
			self.fitline[m].connect()
		
		self.fig.canvas.draw()

	def modify(self,event):
		'''
		Usage:
		To update the projected region based on updated user's input gatefile
		and axes.
		Press 'l' to reload the gatefile and update the axes.
		Follow the instruction from terminal window.
		'''		
		if event.key != 'l':
			return

		self.gatefile=open(gatefile)
		gateline=self.gatefile.readlines()
		self.gate=[]
		for line in gateline:
			if line.find('#') == -1:
				self.gate.append(list(line.split()))
			else:
				self.nuclei=list(line.split())			
		
		self.gatenum=len(self.gate)
		self.gatefile.close()
		print("Enter xlow and xup:\n")
		xlow,xup=map(int,input().split())

		#Transform xlow and xup index in terms of xrebin and yrebin coordinates:
		xlowtrf = np.argmin(np.abs(self.x - xlow))
		xuptrf = np.argmin(np.abs(self.x - xup))

		print("processing...\n")
		print("Which Axes to Update:\n")
		print("Enter axlow and axup:\n")
		axlow,axup=map(int,input().split())


		for i in range(axlow,axup+1,1):
			project=self.project(int(self.gate[i][0]),int(self.gate[i][1]),int(self.gate[i][2]),int(self.gate[i][3]),int(self.gate[i][4]),int(self.gate[i][5]),self.ytrp)
			self.projectrebin[i]=r.rebin(self.xdim,project,self.rebinfactor).rebin()

		for j in range(axlow,axup+1,1):
			if self.gatenum > 1:
				self.ax[j].clear()
				self.ax[j].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				self.ax[j].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				self.ax[j].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.ax[j].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
			else:
				self.ax.clear()
				self.ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				self.ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				self.ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				self.ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
		
		colorlist=['r','g','b','c','k','m','r','g','b','c','k']
		for k in range(axlow,axup+1,1):
			if self.gatenum > 1 :
				lineplotind,=self.ax[k].plot(self.x[xlowtrf:xuptrf],self.projectrebin[k][xlowtrf:xuptrf],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
				self.lineplot[k]=lineplotind
				self.ax[k].legend()
				self.ax[k].set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				self.ax[k].set_ylim(bottom=0)
				self.ax[k].set_xlim(xlow,xup)
			else:
				lineplotind,=self.ax.plot(self.x[xlowtrf:xuptrf],self.projectrebin[k][xlowtrf:xuptrf],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
				self.lineplot[0]=lineplotind
				self.ax.legend()
				self.ax.set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				self.ax.set_ylim(bottom=0)
				self.ax.set_xlim(xlow,xup)
			
		if self.gatenum > 1:
			self.ax[self.gatenum-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
		else:
			self.ax.set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
		
		for m in range(axlow,axup+1,1):
			self.fitline[m]=f.fitgui(self.lineplot[m])
			self.fitline[m].connect()

		self.fig.canvas.draw()
		print("Finished Processing.\n")

	def savefigure(self,event):
		'''
		Usage:
		to save the current figure.	
		'''
		if event.key != 'x':
			return
		plt.get_current_fig_manager()
		plt.savefig(fname=self.nuclei[1]+'.pdf',dpi='figure')


if __name__ == '__main__':
	import sys
	infile=sys.argv[1]
	gatefile=sys.argv[2]
	xlow=int(sys.argv[3])
	xup=int(sys.argv[4])
	rebinfactor=int(sys.argv[5])
	
	gammagate=gggateplotter(infile,5000,5000,0,4999,0,4999,gatefile)
	gammagate.plot(xlow,xup,rebinfactor)
