#!/usr/bin/env python3


import numpy as np
import math

class rebin:
	'''
	Rebin class. Useful for rebinning parsed 1D array.
	07/22 2024
	'''

	def __init__(self,ydim,y,binfactor):
		'''
		Instantiation:
		ydim: the length of parsed 1D array
		y: 1D array
		binfactor: rebinfactor
		'''
		self.ydim=ydim
		self.y=y
		self.binfactor=binfactor

	def rebin(self):
		'''
		Return the bin values, i.e the counts for each bin.
		'''
		dim = int(self.ydim)
		dimb = math.ceil(dim/self.binfactor) #transformed number of dimension after rebinning
		yrebin=np.zeros(int(self.ydim),dtype=np.int32)
		yrebinb=np.zeros(dimb,dtype=np.int32)
		n = 0
		for i in range(0,int(self.ydim),int(self.binfactor)):
			n = n + 1
			for j in range(0,int(self.binfactor),1):
				if i <= int(self.ydim)-int(self.binfactor):
					yrebin[i]=yrebin[i]+self.y[i+j]
			for k in range(0,int(self.binfactor),1):
				if i<= int(self.ydim)-int(self.binfactor):
					yrebin[i+k]=yrebin[i]
			yrebinb[n-1]=yrebin[i]

		return yrebinb

	def rebinx(self):
		'''
		Return the bin center coordinates.
		'''
		p = int(self.binfactor)
		dimb = math.ceil(self.ydim/self.binfactor)
		xrebin=np.zeros(dimb,dtype=np.int32)
		ind = 0
		for m in range(0,int(self.ydim),p):
			ind = ind + 1
			xrebin[ind-1] = m + (p-1)/2.


		return xrebin


if __name__ == "__main__":
	import sys
	import rebin as r
	import matplotlib.pyplot as plt
	import fitgui2 as g

	rebinfactor=int(sys.argv[1])
	z2=np.arange(10,dtype=np.int32)
	z=np.zeros(10,dtype=np.int32)
	z[0]=1
	z[1]=2
	z[2]=3
	z[3]=4
	z[4]=5
	z[5]=5
	z[6]=4
	z[7]=3
	z[8]=2
	z[9]=1
	q1=rebin(10,z,rebinfactor)
	q2=q1.rebin()
	q3=q1.rebinx()
	n1=len(q2)
	n2=len(q3)

	for k in range(n1):
		print("x:",q3[k],"y:",q2[k])

	fig,ax=plt.subplots()
	
	t=r.rebin(10,z,rebinfactor).rebin()
	
	line2, = ax.plot(q3,q2,linewidth=0.85,linestyle='steps-mid')
	fitline2= g.fitgui(line2)
	fitline2.connect()
	plt.show()

