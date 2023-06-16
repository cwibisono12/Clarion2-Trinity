#!/usr/bin/env python3

##Clas for 2D-Rebin:
##Usage: Only for Clarion2-Trinity
#C. Wibisono 06/16/2023

import numpy as np

class rebin2d:

	def __init__(self,ytrp,ydim,xdim,binyfactor,binxfactor):
		self.ytrp=ytrp
		self.ydim=ydim
		self.xdim=xdim
		self.binyfactor=binyfactor
		self.binxfactor=binxfactor

	def rebin(self):
		ytrprebin=np.zeros((int(self.ydim),int(self.xdim)),dtype=np.int32)
		for i in range(0,int(self.ydim),int(self.binyfactor)):
			for m in range(0,int(self.xdim),int(self.binxfactor)):
				for j in range(0,int(self.binyfactor),1):
					for l in range(0,int(self.binxfactor),1):
						if i <= int(self.ydim)-int(self.binyfactor):
							if m <= int(self.xdim)-int(self.binxfactor):
								ytrprebin[i][m]=ytrprebin[i][m]+self.ytrp[i+j][m+l]

				for k in range(0,int(self.binyfactor),1):
					for n in range(0,int(self.binxfactor),1):
						if i<= int(self.ydim)-int(self.binyfactor):
							if m <= int(self.xdim)-int(self.binxfactor):
								ytrprebin[i+k][m+n]=ytrprebin[i][m]

		return ytrprebin
