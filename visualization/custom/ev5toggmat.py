#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np

def ev5toggmat(fpr,*,dimy=5000,dimx=5000):
	'''
	Function to Generate ggmatrix from ev5 file
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension
	dimx: x dimension
	'''
	#Prompt:
	ggmatp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp_nd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp_dd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#Non-Prompt:
	ggmatnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatnp_nd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatnp_dd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)

	#tdiff:
	ggmat_tdiff=np.ndarray(shape=(dimy,dimx),dtype=np.int32)

	for i in range(dimy):
		for j in range(dimx):
			ggmatp[i][j]=0
			ggmatp_nd[i][j]=0
			ggmatp_dd[i][j]=0
			ggmatnp[i][j]=0
			ggmatnp_nd[i][j]=0
			ggmatnp_dd[i][j]=0
			ggmat_tdiff[i][j]=0

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if temp[0][0] == 0 and temp[0][1] == 0:
					gmult=temp[0][2]
					for i in temp[1].keys():
						for j in temp[1].keys():
							if i!=j :
								tdiff=temp[1][i][6] - temp[1][j][6] + 2000
								if (tdiff > 0 and tdiff <5000) and (temp[1][i][2] >=0 and temp[1][i][2] <5000)	and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
									if temp[1][i][2] < temp[1][j][2]:
										ggmat_tdiff[temp[1][i][2]][tdiff] = ggmat_tdiff[temp[1][i][2]][tdiff] + 1	
									else:											
										ggmat_tdiff[temp[1][j][2]][tdiff] = ggmat_tdiff[temp[1][i][2]][tdiff] + 1	



								if tdiff >= 1967 and tdiff <=2032:
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
										ggmatp[temp[1][i][2]][temp[1][j][2]]=ggmatp[temp[1][i][2]][temp[1][j][2]] + 1	
									if (temp[1][i][3] >= 0 and temp[1][i][3] <5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
										ggmatp_nd[temp[1][i][3]][temp[1][j][2]]=ggmatp_nd[temp[1][i][3]][temp[1][j][2]] + 1
									if (temp[1][i][3] >= 0 and temp[1][i][3] < 5000) and (temp[1][j][3] >=0 and temp[1][j][3] <5000):
										ggmatp_dd[temp[1][i][3]][temp[1][j][3]]=ggmatp_dd[temp[1][i][3]][temp[1][j][3]] + 1


								if (tdiff >= 1840 and tdiff <=1872) or (tdiff >= 2130 and tdiff <-2162) :
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
										ggmatnp[temp[1][i][2]][temp[1][j][2]]=ggmatnp[temp[1][i][2]][temp[1][j][2]] + 1	
									if (temp[1][i][3] >= 0 and temp[1][i][3] <5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
										ggmatnp_nd[temp[1][i][3]][temp[1][j][2]]=ggmatnp_nd[temp[1][i][3]][temp[1][j][2]] + 1
									if (temp[1][i][3] >= 0 and temp[1][i][3] < 5000) and (temp[1][j][3] >=0 and temp[1][j][3] <5000):
										ggmatnp_dd[temp[1][i][3]][temp[1][j][3]]=ggmatnp_dd[temp[1][i][3]][temp[1][j][3]] + 1



	matwrite("ggmat_tdiff.spn2",dimy=dimy,dimx=dimx,arr=ggmat_tdiff,overwrite=1)
	matwrite("ggmatp.spn2",dimy=dimy,dimx=dimx,arr=ggmatp,overwrite=1)
	matwrite("ggmatp_nd.spn2",dimy=dimy,dimx=dimx,arr=ggmatp_nd,overwrite=1)
	matwrite("ggmatp_dd.spn2",dimy=dimy,dimx=dimx,arr=ggmatp_dd,overwrite=1)
	matwrite("ggmatnp.spn2",dimy=dimy,dimx=dimx,arr=ggmatnp,overwrite=1)
	matwrite("ggmatnp_nd.spn2",dimy=dimy,dimx=dimx,arr=ggmatnp_nd,overwrite=1)
	matwrite("ggmatnp_dd.spn2",dimy=dimy,dimx=dimx,arr=ggmatnp_dd,overwrite=1)


if __name__ == "__main__":
	import sys
	fileev5=sys.argv[1]
	ev5toggmat(fileev5,dimy=5000,dimx=5000)
