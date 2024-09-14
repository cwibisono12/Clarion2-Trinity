#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread 

def ev5toadg(fpr,xlow,xup,xlowbkg,xupbkg,*,dimy=16,dimx=8192,overwrite=0):
	'''
	Function to Generate Angular Distribution Spectra from ev5 file gated on particle and gamma.
	Parameter(s):
	fpr: ev5 file pointer object
	xlow: int lower limit of gamma gate
	xup: int upper limit of gamma gate
	xlowbkg: int lower limit of background region
	xup: int upper limit of background region
	dimy: y dimension
	dimx: x dimension
	overwrite: 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	p1_angularp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p1_angularp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularp_dopp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularp_dopp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#Non-Prompt:
	p1_angularnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularnp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularnp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p1_angularnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularnp_dopp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_angularnp_dopp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	for i in range(dimy):
		for j in range(dimx):
			p1_angularp[i][j]=0
			p1_angularp_p[i][j]=0
			p1_angularp_np[i][j]=0
			
			p1_angularp_dopp[i][j]=0
			p1_angularp_dopp_p[i][j]=0
			p1_angularp_dopp_np[i][j]=0
			
			p1_angularnp[i][j]=0
			p1_angularnp_p[i][j]=0
			p1_angularnp_np[i][j]=0
			
			p1_angularnp_dopp[i][j]=0
			p1_angularnp_dopp_p[i][j]=0
			p1_angularnp_dopp_np[i][j]=0

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if temp[0][0] == 1 and temp[0][1] == 0: #1p-gated
					for j in temp[1].keys():
						'''Prompt Region'''
						if ((temp[1][j][2] >= xlow and temp[1][j][2] <= xup) and temp[1][j][4] == 1):
							for i in temp[1].keys():
								if (i!=j):
									if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][4] == 1):
										p1_angularp_p[int(i)][temp[1][i][2]]=p1_angularp_p[int(i)][temp[1][i][2]]+1
									if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][4] == 1):
										p1_angularp_dopp_p[int(i)][temp[1][i][3]]=p1_angularp_dopp_p[int(i)][temp[1][i][3]]+1
							
					
									if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][5] == 1):
										p1_angularnp_p[int(i)][temp[1][i][2]]=p1_angularnp_p[int(i)][temp[1][i][2]]+1
									if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][5] == 1):
										p1_angularnp_dopp_p[int(i)][temp[1][i][3]]=p1_angularnp_dopp_p[int(i)][temp[1][i][3]]+1
						'''Background Region'''
						if ((temp[1][j][2] >= xlowbkg and temp[1][j][2] <= xupbkg) and temp[1][j][4] == 1):
							for i in temp[1].keys():
								if (i!=j):
									if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][4] == 1):
										p1_angularp_np[int(i)][temp[1][i][2]]=p1_angularp_np[int(i)][temp[1][i][2]]+1
									if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][4] == 1):
										p1_angularp_dopp_np[int(i)][temp[1][i][3]]=p1_angularp_dopp_np[int(i)][temp[1][i][3]]+1
							
					
									if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][5] == 1):
										p1_angularnp_np[int(i)][temp[1][i][2]]=p1_angularnp_np[int(i)][temp[1][i][2]]+1
									if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][5] == 1):
										p1_angularnp_dopp_np[int(i)][temp[1][i][3]]=p1_angularnp_dopp_np[int(i)][temp[1][i][3]]+1
		
			
	for k1 in range(dimy):
		for k2 in range(dimx):
			p1_angularp[k1][k2] = p1_angularp_p[k1][k2] - p1_angularp_np[k1][k2]
			p1_angularp_dopp[k1][k2] = p1_angularp_dopp_p[k1][k2] - p1_angularp_dopp_np[k1][k2]
			p1_angularnp[k1][k2] = p1_angularnp_p[k1][k2] - p1_angularnp_np[k1][k2]
			p1_angularnp_dopp[k1][k2] = p1_angularnp_dopp_p[k1][k2] - p1_angularnp_dopp_np[k1][k2]

	t1 = Thread(target=matwrite,args=("p1_angularp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_angularp,'overwrite':overwrite})
	t2 = Thread(target=matwrite,args=("p1_angularp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_angularp_dopp,'overwrite':overwrite})
	
	t1b = Thread(target=matwrite,args=("p1_angularnp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_angularnp,'overwrite':overwrite})
	t2b = Thread(target=matwrite,args=("p1_angularnp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_angularnp_dopp,'overwrite':overwrite})
	
	t1.start()
	t2.start()

	t1b.start()
	t2b.start()

	t1.join()
	t2.join()

	t1b.join()
	t2b.join()


if __name__ == "__main__":
	import sys
	fileev5 = sys.argv[1]
	overwrite = int(sys.argv[2])
	xlow = int(sys.argv[3])
	xup = int(sys.argv[4])
	xlowbkg = int(sys.argv[5])
	xupbkg = int(sys.argv[6])

	ev5toadg(fileev5,xlow,xup,xlowbkg,xupbkg,dimy=16,dimx=8192,overwrite=overwrite)
