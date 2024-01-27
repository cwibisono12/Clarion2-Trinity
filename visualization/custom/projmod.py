#!/usr/bin/env python3

###Class for Projection and Histogram Reader
##2nd Version-- 06/11/2023
###Usage: Only for Clarion2-Trinity
#C. Wibisono 10/19/2022


import numpy as np

class clarion:
	
	def __init__(self,infile,ydim,xdim,xlow,xup,ylow,yup):
		self.infile=infile
		self.ydim=ydim
		self.xdim=xdim
		self.xlow=xlow
		self.xup=xup
		self.ylow=ylow
		self.yup=yup
	
	def readparse(self):
		y=np.fromfile(self.infile,dtype=np.int32,sep="",count=-1)
		ytrp=np.reshape(y,(int(self.ydim),int(self.xdim)))
		return ytrp


	def project(self,proj,background,widthleft,widthright,backleft,backright,ytrp):
		xproj=np.zeros(int(self.xdim),dtype=np.int32)
		yproj=np.zeros(int(self.ydim),dtype=np.int32)
		xprojb=np.zeros(int(self.xdim),dtype=np.int32)
		yprojb=np.zeros(int(self.ydim),dtype=np.int32)
		xprojbs=np.zeros(int(self.xdim),dtype=np.int32)
		yprojbs=np.zeros(int(self.ydim),dtype=np.int32)
		#ytrp=self.readparse()
		widthgate=int(widthright)-int(widthleft)
		widthbackg=int(backright)-int(backleft)	
		if int(proj) == 1:
			for i in range(0,int(self.xdim),1):
				for j in range(int(widthleft),int(widthright)+1,1):
					xproj[i]=xproj[i]+ytrp[j,i]
			
			if int(background)==0:
				return xproj
			if int(background)==1:
				for i in range(0,int(self.xdim),1):
					for j in range(int(backleft),int(backright)+1,1):
						xprojb[i]=xprojb[i]+ytrp[j,i]
					xprojbs[i]=(widthgate/(widthgate+widthbackg))*xproj[i]-(widthbackg/(widthgate+widthbackg))*xprojb[i]
				return xprojbs
		
	
		if int(proj) == 2:
			for j in range(0,int(self.ydim),1):
				for i in range(int(widthleft),int(widthright)+1,1):
					yproj[j]=yproj[j]+ytrp[j,i]
			if int(background)==0:
				return yproj
			if int(background)==1:
				for j in range(0,int(self.ydim),1):
					for i in range(int(backleft),int(backright)+1,1):
						yprojb[j]=yprojb[j]+ytrp[j,i]
		
					yprojbs[j]=(widthgate/(widthgate+widthbackg))*yproj[j]-(widthbackg/(widthgate+widthbackg))*yprojb[j]
				return yprojbs
	

	
