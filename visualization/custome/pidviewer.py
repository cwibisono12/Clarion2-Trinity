#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import projmod as p
from matplotlib.colors import LogNorm
import cut2D as cut

class pidviewer:
	'''
	PID Viewer for Clarion2-Trinity
	
	Usage:
	To view the entire PID's of Clarion2-Trinity 
	This Class is set to be system's dependent. In other words, the use for another detector's systems might 
	require several adjustments for the class' methods.

	C. Wibisono
	10/23 '23
	'''
	plt.rcParams['font.family']='serif'
	plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

	def __init__(self,pidlist,dimy,dimx,*,bangate=None):
		'''	
		Constructor Description:
		1). pidlist--->list of pid file
		2). dimy---> y dimension of a pid
		3). dimx---> x dimension of a pid
		4). bangate(optional) ----> banana cuts file if we have those already.
		'''
		pidline=open(pidlist).readlines()
		self.pid=[]
		b=[]
		row=0
		for line in pidline:
			if line.find('#') == -1:
				self.pid.append(line.split())
				b.append(self.pid[row][1])
				row=row+1
			else:
				pass

		self.ring=list(set(b))
		self.idring=[]
		del b
		#Count the number of id for each ring:
		for i in range(len(self.ring)):
			count=0
			for j in range(len(self.pid)):
				if self.ring[i] == self.pid[j][1]:
					count=count+1
			self.idring.append(count)	

		self.dimy=dimy
		self.dimx=dimx

		if bangate != None:
			self.bangate=bangate

	def plot(self,xlow,xup,ylow,yup):
		'''
		Method's Description:
		To plot PID's all at once based on the given pidlist files.
		Arguments:
		1.xlow--->lower limit of x axis
		2.xup--->upper limit of x axis
		3.ylow-->lower limit of y axis
		4.yup--->upper limit of y axis
		'''			
		for i in range(len(self.ring)):
			if (self.idring[i] % 4 ==0) and (self.idring[i] % 2 ==0):
				if self.ring[i] == '1':
					fig1, ax1=plt.subplots(int(self.idring[i]/4),4)
					fig1.suptitle('Ring 1\nGAGG:Ce-Si-Si')
				if self.ring[i] == '4':
					fig4, ax4=plt.subplots(int(self.idring[i]/4),4)	
					fig4.suptitle('Ring 4\nGAGG:Ce-Si-Si')
				if self.ring[i] == '5':
					fig5, ax5=plt.subplots(int(self.idring[i]/4),4)
					fig5.suptitle('Ring 5\nGAGG:Ce-Si-Si')
			if (self.idring[i] % 2 ==0) and (self.idring[i] % 4 !=0):
				if self.ring[i] == '2':
					fig2, ax2=plt.subplots(2,int(self.idring[i]/2))
					fig2.suptitle('Ring 2\nGAGG:Ce-Si-Si')
				if self.ring[i] == '3':
					fig3, ax3=plt.subplots(2,int(self.idring[i]/2))
					fig3.suptitle('Ring 3\nGAGG:Ce-Si-Si')
		
		for i in range(len(self.ring)):
			if self.ring[i] == '1':
				for m1 in range(int(self.idring[i]/4)):
					for m2 in range(4):
						ax1[m1,m2].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.7)
						ax1[m1,m2].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.7)
						ax1[m1,m2].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
						ax1[m1,m2].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
				
			if self.ring[i] == '2':
				for m1 in range(2):
					for m2 in range(int(self.idring[i]/2)):
						ax2[m1,m2].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.7)
						ax2[m1,m2].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.7)
						ax2[m1,m2].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
						ax2[m1,m2].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

			if self.ring[i] == '3':
				for m1 in range(2):
					for m2 in range(int(self.idring[i]/2)):
						ax3[m1,m2].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.7)
						ax3[m1,m2].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.7)
						ax3[m1,m2].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
						ax3[m1,m2].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

			if self.ring[i] == '4':
				for m1 in range(int(self.idring[i]/4)):
					for m2 in range(4):
						ax4[m1,m2].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.7)
						ax4[m1,m2].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.7)
						ax4[m1,m2].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
						ax4[m1,m2].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

			if self.ring[i] == '5':
				for m1 in range(int(self.idring[i]/4)):
					for m2 in range(4):
						ax5[m1,m2].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.7)
						ax5[m1,m2].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.7)
						ax5[m1,m2].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
						ax5[m1,m2].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

		
		
				
			
		for m in range(len(self.pid)):
			obj=p.clarion(self.pid[m][0],self.dimy,self.dimx,0,self.dimx-1,0,self.dimy-1)
			pidobj=obj.readparse()
			#objrebin=r.rebin2d(pidobj,self.dimy,self.dimx,self.rebiny,self.rebinx).rebin()
			if self.pid[m][1] == '1':
				if int(self.pid[m][2])-1>=0 and int(self.pid[m][2])-1<=3:
					ax1[0,int(self.pid[m][2])-1].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax1[0,int(self.pid[m][2])-1].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=4 and int(self.pid[m][2])-1<=7:
					ax1[1,int(self.pid[m][2])-1-4].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax1[1,int(self.pid[m][2])-1-4].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])

			if self.pid[m][1] == '2':
				if int(self.pid[m][2])-1>=0 and int(self.pid[m][2])-1<=4:
					ax2[0,int(self.pid[m][2])-1].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax2[0,int(self.pid[m][2])-1].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				else:
					ax2[1,int(self.pid[m][2])-1-5].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax2[1,int(self.pid[m][2])-1-5].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])

			if self.pid[m][1] == '3':
				if int(self.pid[m][2])-1>=0 and int(self.pid[m][2])-1<=6:
					ax3[0,int(self.pid[m][2])-1].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax3[0,int(self.pid[m][2])-1].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=7 and int(self.pid[m][2])-1<=13:
					ax3[1,int(self.pid[m][2])-1-7].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax3[1,int(self.pid[m][2])-1-7].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])

			if self.pid[m][1] == '4':
				if int(self.pid[m][2])-1>=0 and int(self.pid[m][2])-1<=3:
					ax4[0,int(self.pid[m][2])-1].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax4[0,int(self.pid[m][2])-1].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=4 and int(self.pid[m][2])-1<=7:
					ax4[1,int(self.pid[m][2])-1-4].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax4[1,int(self.pid[m][2])-1-4].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=8 and int(self.pid[m][2])-1<=11:
					ax4[2,int(self.pid[m][2])-1-8].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax4[2,int(self.pid[m][2])-1-8].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=12 and int(self.pid[m][2])-1<=16:
					ax4[3,int(self.pid[m][2])-1-12].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax4[3,int(self.pid[m][2])-1-12].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])

			if self.pid[m][1] == '5':
				if int(self.pid[m][2])-1>=0 and int(self.pid[m][2])-1<=3:
					ax5[0,int(self.pid[m][2])-1].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax5[0,int(self.pid[m][2])-1].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=4 and int(self.pid[m][2])-1<=7:
					ax5[1,int(self.pid[m][2])-1-4].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax5[1,int(self.pid[m][2])-1-4].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=8 and int(self.pid[m][2])-1<=11:
					ax5[2,int(self.pid[m][2])-1-8].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax5[2,int(self.pid[m][2])-1-8].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
				if int(self.pid[m][2])-1>=12 and int(self.pid[m][2])-1<=16:
					ax5[3,int(self.pid[m][2])-1-12].imshow(pidobj,cmap='gist_ncar',origin='lower',extent=[0,4095,0,4095],aspect='auto',norm=LogNorm())
					ax5[3,int(self.pid[m][2])-1-12].set_title(label=self.pid[m][1]+'-'+self.pid[m][2])
		

		for i in range(len(self.ring)):
			if self.ring[i] == '1':
				p1=[]
				for m1 in range(int(self.idring[i]/4)):
					for m2 in range(4):
						p1.append(cut.cut2D(ax1[m1,m2],bangate=self.bangate))
						ax1[m1,m2].set_xlim(xlow,xup)
						ax1[m1,m2].set_ylim(ylow,yup)
				for l in range(len(p1)):
					p1[l].connect()
					p1[l].cutfile(l)
	
			if self.ring[i] == '2':
				p2=[]
				for m1 in range(2):
					for m2 in range(int(self.idring[i]/2)):
						p2.append(cut.cut2D(ax2[m1,m2],bangate=self.bangate))
						ax2[m1,m2].set_xlim(xlow,xup)
						ax2[m1,m2].set_ylim(ylow,yup)
				for l in range(len(p2)):
					p2[l].connect()
					p2[l].cutfile(l+1)

			if self.ring[i] == '3':
				p3=[]
				for m1 in range(2):
					for m2 in range(int(self.idring[i]/2)):
						p3.append(cut.cut2D(ax3[m1,m2],bangate=self.bangate))
						ax3[m1,m2].set_xlim(xlow,xup)
						ax3[m1,m2].set_ylim(ylow,yup)
				for l in range(len(p3)):
					p3[l].connect()
					p3[l].cutfile(l)

			if self.ring[i] == '4':
				p4=[]
				for m1 in range(int(self.idring[i]/4)):
					for m2 in range(4):
						p4.append(cut.cut2D(ax4[m1,m2],bangate=self.bangate))
						ax4[m1,m2].set_xlim(xlow,xup)
						ax4[m1,m2].set_ylim(ylow,yup)
				for l in range(len(p4)):
					p4[l].connect()
					p4[l].cutfile(l+11)

			if self.ring[i] == '5':
				p5=[]
				for m1 in range(int(self.idring[i]/4)):
					for m2 in range(4):
						p5.append(cut.cut2D(ax5[m1,m2],bangate=self.bangate))
						ax5[m1,m2].set_xlim(xlow,xup)
						ax5[m1,m2].set_ylim(ylow,yup)

				for l in range(len(p5)):
					p5[l].connect()
					p5[l].cutfile(l)

		plt.show()


if __name__ == "__main__":
	import sys
	pidlist=sys.argv[1]
	dimy=int(sys.argv[2])
	dimx=int(sys.argv[3])
	banfile=sys.argv[4]
	gagg=pidviewer(pidlist,dimy,dimx,bangate=banfile)
	gagg.plot(1400,3000,0,3000)
