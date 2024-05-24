#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np


def ev5toad(fpr,*,dimy=16,dimx=8192,overwrite=0):
	'''
	Function to Generate Angular Distribution Spectra from ev5 file
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension
	dimx: x dimension
	overwrite: 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	p2_angularp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_angularp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	#p3_angularp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p2_angularp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_angularp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	#p3_angularp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#Non-Prompt:
	p2_angularnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_angularnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	#p3_angularnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p2_angularnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_angularnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	#p3_angularnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	for i in range(dimy):
		for j in range(dimx):
			p2_angularp[i][j]=0
			pa_angularp[i][j]=0
			#p3_angularp[i][j]=0
			
			p2_angularp_dopp[i][j]=0
			pa_angularp_dopp[i][j]=0
			#p3_angularp_dopp[i][j]=0
			
			p2_angularnp[i][j]=0
			pa_angularnp[i][j]=0
			#p3_angularnp[i][j]=0
			
			p2_angularnp_dopp[i][j]=0
			pa_angularnp_dopp[i][j]=0
			#p3_angularnp_dopp[i][j]=0

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if temp[0][0] == 2 and temp[0][1] == 0: #2p-gated
					for i in temp[1].keys():
						if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][4] == 1):
							p2_angularp[int(i)][temp[1][i][2]]=p2_angularp[int(i)][temp[1][i][2]]+1
						if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][4] == 1):
							p2_angularp_dopp[int(i)][temp[1][i][3]]=p2_angularp_dopp[int(i)][temp[1][i][3]]+1
						if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][5] == 1):
							p2_angularnp[int(i)][temp[1][i][2]]=p2_angularnp[int(i)][temp[1][i][2]]+1
						if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][5] == 1):
							p2_angularnp_dopp[int(i)][temp[1][i][3]]=p2_angularnp_dopp[int(i)][temp[1][i][3]]+1

				if temp[0][0] == 1 and temp[0][1] == 1: #1p-1agated
					for i in temp[1].keys():
						if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][4] == 1):
							pa_angularp[int(i)][temp[1][i][2]]=pa_angularp[int(i)][temp[1][i][2]]+1
						if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][4] == 1):
							pa_angularp_dopp[int(i)][temp[1][i][3]]=pa_angularp_dopp[int(i)][temp[1][i][3]]+1
						if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][5] == 1):
							pa_angularnp[int(i)][temp[1][i][2]]=pa_angularnp[int(i)][temp[1][i][2]]+1
						if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][5] == 1):
							pa_angularnp_dopp[int(i)][temp[1][i][3]]=pa_angularnp_dopp[int(i)][temp[1][i][3]]+1

			

	matwrite("p2_angularp.spn2",dimy=dimy,dimx=dimx,arr=p2_angularp,overwrite=overwrite)	
	matwrite("p2_angularp_dopp.spn2",dimy=dimy,dimx=dimx,arr=p2_angularp_dopp,overwrite=overwrite)	
	matwrite("pa_angularp.spn2",dimy=dimy,dimx=dimx,arr=pa_angularp,overwrite=overwrite)	
	matwrite("pa_angularp_dopp.spn2",dimy=dimy,dimx=dimx,arr=pa_angularp_dopp,overwrite=overwrite)	
	
	matwrite("p2_angularnp.spn2",dimy=dimy,dimx=dimx,arr=p2_angularnp,overwrite=overwrite)	
	matwrite("p2_angularnp_dopp.spn2",dimy=dimy,dimx=dimx,arr=p2_angularnp_dopp,overwrite=overwrite)	
	matwrite("pa_angularnp.spn2",dimy=dimy,dimx=dimx,arr=pa_angularnp,overwrite=overwrite)	
	matwrite("pa_angularnp_dopp.spn2",dimy=dimy,dimx=dimx,arr=pa_angularnp_dopp,overwrite=overwrite)	


if __name__ == "__main__":
	import sys
	fileev5 = sys.argv[1]
	overwrite = int(sys.argv[2])

	ev5toad(fileev5,dimy=16,dimx=8192,overwrite=overwrite)
