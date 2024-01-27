#!/usr/bin/env python3

##Clas for Rebin:
##Usage: Only for Clarion2-Trinity
#C. Wibisono 02/25/2023

import numpy as np

class rebin:

	def __init__(self,ydim,y,binfactor):
		self.ydim=ydim
		self.y=y
		self.binfactor=binfactor

	def rebin(self):
		yrebin=np.zeros(int(self.ydim),dtype=np.int32)
		for i in range(0,int(self.ydim),int(self.binfactor)):
			for j in range(0,int(self.binfactor),1):
				if i <= int(self.ydim)-int(self.binfactor):
					yrebin[i]=yrebin[i]+self.y[i+j]
			for k in range(0,int(self.binfactor),1):
				if i<= int(self.ydim)-int(self.binfactor):
					yrebin[i+k]=yrebin[i]

		return yrebin
