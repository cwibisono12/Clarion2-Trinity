#!/usr/bin/env python3

import numpy as np
import projmod as p
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import fitgui as f
import rebin as r

class genmatplotter:
	'''
	
	Class Description:
	Custome Class to visualize projections from 2D-matrices given prompt and non-prompt matrices and list of gates.
	Intended for Angular Distribution Fit, Polarization Fit, and many more related with gamma spectroscopy.
	
	The use case is not necessarily limited to above cases. 
	
	C. Wibisono
	10/17 '23
	
	'''	
	plt.rcParams['font.family']='serif'
	plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

	def __init__(self,matfilelist,gatefile):
		'''
		Constructor Description:
		1).matfilelist ---> a file consisting of lists of prompt and non-prompt matrices followed by
					dimensional information
				    see:Clarion2-Trinity/param/matfile.txt as an example
		2).gatefile ----> a file consisting of list of gates
				    see:Clarion2-Trinity/param/gatefilead.txt as an example
		'''
		matfile=open(matfilelist)
		matlinelist=matfile.readlines()
		self.matlist=[]
		for line in matlinelist:
			if line.find('#') == -1:
				self.matlist.append(list(line.split()))
			else:
				pass

		self.matnum=len(self.matlist)
		gatefiles=open(gatefile)
		gateline=gatefiles.readlines()
		self.gate=[]
		for line in gateline:
			if line.find('#') == -1:
				self.gate.append(list(line.split()))
			else:
				self.nuclei=list(line.split())
		
		self.gatenum=len(self.gate)
		#print(int(self.matlist[0][3]))
		self.x=np.arange(0,int(self.matlist[0][3]),1)
		self.project=[]
		self.projectrebin=[]
		self.flagad=0
		matfile.close()

	def process(self,rebinfactor):
		'''
		Method's Description:
		To project a matrix from the list of gates
		'''
		pmat=[] #List of prompt matrices
		npmat=[] #List of non-prompt matrices
		ptrp=[] #List of parsed prompt-matrices
		nptrp=[] #List of parsed non-prompt-matrices

		for i in range(self.matnum):
			if self.matlist[i][0] == 'p':
				zpmat=p.clarion(self.matlist[i][1],int(self.matlist[i][2]),int(self.matlist[i][3]),int(self.matlist[i][4]),int(self.matlist[i][5]),int(self.matlist[i][6]),int(self.matlist[i][7]))
				zptrp=zpmat.readparse()
				pmat.append(zpmat)
				ptrp.append(zptrp)
			if self.matlist[i][0] == 'np':
				znpmat=p.clarion(self.matlist[i][1],int(self.matlist[i][2]),int(self.matlist[i][3]),int(self.matlist[i][4]),int(self.matlist[i][5]),int(self.matlist[i][6]),int(self.matlist[i][7]))
				znptrp=znpmat.readparse()
				npmat.append(znpmat)
				nptrp.append(znptrp)
		
		#Projected Parsed Matrices:
		for k in range(self.gatenum):
			for l in range(len(pmat)):
				pproj=pmat[l].project(int(self.gate[k][0]),int(self.gate[k][1]),int(self.gate[k][2]),int(self.gate[k][3]),int(self.gate[k][4]),int(self.gate[k][5]),ptrp[l])
				npproj=npmat[l].project(int(self.gate[k][0]),int(self.gate[k][1]),int(self.gate[k][2]),int(self.gate[k][3]),int(self.gate[k][4]),int(self.gate[k][5]),nptrp[l])
				projclean=pproj-npproj
				self.project.append(projclean)
				#self.projectrebin.append(r.rebin(int(self.matlist[0][3]),self.project[k],rebinfactor).rebin())
		for n in range(len(self.project)):
			self.projectrebin.append(r.rebin(int(self.matlist[0][3]),self.project[n],rebinfactor).rebin())

		self.pmatlen=len(pmat) #number of prompt matrices

	def ad(self,rebinfactor):
		'''
		Method's Description:
		To reprocess lists of gates and stack the spectra according to the group of Clarion ring's angles
		'''
		self.addforw=np.zeros(int(self.matlist[0][3]))
		self.add90=np.zeros(int(self.matlist[0][3]))
		self.addback1=np.zeros(int(self.matlist[0][3]))
		self.addback2=np.zeros(int(self.matlist[0][3]))

		for i in range(self.gatenum):
			if self.gate[i][7] == 'forward':
				self.addforw=self.addforw+self.project[i]
			if self.gate[i][7] == '90':
				self.add90=self.add90+self.project[i]
			if self.gate[i][7] == 'backward1':
				self.addback1=self.addback1+self.project[i]
			if self.gate[i][7] == 'backward2':
				self.addback2=self.addback2+self.project[i]
			
		del self.project
		del self.projectrebin
		self.project=[self.addforw,self.add90,self.addback1,self.addback2]
		self.projectrebin=[]
		for k in range(len(self.project)):
			self.projectrebin.append(r.rebin(int(self.matlist[0][3]),self.project[k],rebinfactor).rebin())
			
		self.flagad=1
	
		
	def plot(self,xlow,xup,title):
		'''
		Method's Description:
		To visualize the projection
		'''
		if self.flagad == 0:
			iterator=len(self.project)
		if self.flagad == 1:
			iterator=len(self.project)

		fig,ax=plt.subplots(iterator,1)
		fig.suptitle(title)

		labelad=['48.24','90','131.75','150']
		labelnonad=[]

		for k1 in range(self.gatenum):
			labk1=self.gate[k1][6]
			for k2 in range(self.pmatlen):
				labelnonad.append(labk1)

		for m in range(iterator):
			if iterator > 1:
				ax[m].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				ax[m].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				ax[m].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				ax[m].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
			else:
				ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
				ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
				ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	

		lineplot=[]
		colorlist=['r','g','b','c','k','m','r','g','b','c','k']
		for n in range(iterator):
			if self.flagad == 0:
				if iterator > 1: #temp
					lineplotind,=ax[n].plot(self.x[xlow:xup],self.projectrebin[n][xlow:xup],linewidth=0.85,ls='steps-mid',color=colorlist[n],label=labelnonad[n])
				else:
					lineplotind,=ax.plot(self.x[xlow:xup],self.projectrebin[n][xlow:xup],linewidth=0.85,ls='steps-mid',color=colorlist[n],label=labelnonad[n])

			if self.flagad == 1:
				lineplotind,=ax[n].plot(self.x[xlow:xup],self.projectrebin[n][xlow:xup],linewidth=0.85,ls='steps-mid',color=colorlist[n],label=labelad[n])

			lineplot.append(lineplotind)
			if iterator > 1:
				ax[n].legend()
				ax[n].set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				ax[n].set_ylim(bottom=0)
				ax[n].set_xlim(xlow,xup)
			else:
				ax.legend()
				ax.set_ylabel('counts/'+str(rebinfactor)+' '+'keV',style='normal',fontweight='bold')
				ax.set_ylim(bottom=0)
				ax.set_xlim(xlow,xup)

		if self.flagad == 0:
			if iterator > 1:
				#ax[self.gatenum-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
				ax[iterator-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')
			else:
				ax.set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')

		if self.flagad == 1:
			ax[iterator-1].set_xlabel('E$_{\gamma}$ (keV)',style='normal',fontweight='bold')

		fitline=[]
		for s in range(iterator):
			fitline.append(f.fitgui(lineplot[s]))
			fitline[s].connect()

		plt.show()
	

if __name__ == "__main__":
	import sys
	matfile=sys.argv[1]
	gatefile=sys.argv[2]
	xlow=int(sys.argv[3])
	xup=int(sys.argv[4])
	rebinfactor=int(sys.argv[5])
	option=int(sys.argv[6])
	title="".join(sys.argv[7])

	adfit=genmatplotter(matfile,gatefile)
	adfit.process(rebinfactor)
	if option == 1:
		adfit.ad(rebinfactor)
	adfit.plot(xlow,xup,title)
