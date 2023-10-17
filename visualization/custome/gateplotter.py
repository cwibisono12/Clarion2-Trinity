#!/usr/bin/env python3

import projmod as p
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as tck
import fitgui as f
import rebin as r

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
				self.gate.append(list(map(int,line.split())))
			else:
				self.nuclei=list(line.split())			


		self.gatenum=len(self.gate)
		self.x=np.arange(0,xdim,1)

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
		projectrebin=[]
		
		for i in range(self.gatenum):
			project.append(self.project(self.gate[i][0],self.gate[i][1],self.gate[i][2],self.gate[i][3],self.gate[i][4],self.gate[i][5],self.ytrp))
			projectrebin.append(r.rebin(self.xdim,project[i],rebinfactor).rebin())

		fig,ax=plt.subplots(self.gatenum,1)
		fig.suptitle(str(self.nuclei[1])+' '+'Gamma Spectra'+'\n'+'Clarion2-Trinity')
		for j in range(self.gatenum):
			ax[j].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
			ax[j].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
			ax[j].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
			ax[j].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
					
		lineplot=[]
		colorlist=['r','g','b','c','k','m','r','g','b','c','k']
		for k in range(0,self.gatenum,1):
			lineplotind,=ax[k].plot(self.x[xlow:xup],projectrebin[k][xlow:xup],linewidth=0.85,ls='steps-mid',color=colorlist[k],label=str(self.gate[k][6])+' '+'keV'+' '+'gate')
			lineplot.append(lineplotind)
			#lineplot.append(ax[k].plot(self.x[xlow:xup],projectrebinr[k][xlow:xup],linewidth=0.85,ls='steps-mid',color='r',label=str(self.gate[k][6])))
			ax[k].legend()
			ax[k].set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
			ax[k].set_ylim(bottom=0)
			ax[k].set_xlim(xlow,xup)

			
		ax[self.gatenum-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')

		
		fitline=[]
		for m in range(0,self.gatenum,1):
			fitline.append(f.fitgui(lineplot[m]))
			fitline[m].connect()
		
		plt.show()


if __name__ == '__main__':
	import sys
	infile=sys.argv[1]
	gatefile=sys.argv[2]
	xlow=int(sys.argv[3])
	xup=int(sys.argv[4])
	rebinfactor=int(sys.argv[5])
	
	gammagate=gggateplotter(infile,5000,5000,0,4999,0,4999,gatefile)
	gammagate.plot(xlow,xup,rebinfactor)
