#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread

def ev5toggmatmultqc(fpr,*,dimy=5000,dimx=5000,overwrite=0):
	'''
	Function to Generate ggmatrix from ev5 file
	This version would generate gamma-gamma matrix with one ring angle for one axis and all angles for the other 
	axis.
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension
	dimx: x dimension
	overwrite: overwrite flag 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	ggmatp2_f=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp2_perp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp2_b1=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp2_b2=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	


	for i in range(dimy):
		for j in range(dimx):
			ggmatp2_f[i][j]=0
			ggmatp2_perp[i][j]=0
			ggmatp2_b1[i][j]=0
			ggmatp2_b2[i][j]=0


	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if temp[0][0] == 2 and temp[0][1] == 0: #2p-gated
					gmult=temp[0][2]
					for i in temp[1].keys():
						for j in temp[1].keys():
							if i!=j :
								tdiff=temp[1][i][6] - temp[1][j][6] + 2000
								if tdiff >= 1988 and tdiff <=2012:
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1) and (i == '5' or i == '6' or i == '11'):
										ggmatp2_f[temp[1][i][2]][temp[1][j][2]]=ggmatp2_f[temp[1][i][2]][temp[1][j][2]] + 1	
									if (temp[1][i][2] >= 0 and temp[1][i][2] <5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1) and (i == '3' or i == '4' or i == '9' or i == '10'):
										ggmatp2_perp[temp[1][i][2]][temp[1][j][2]]=ggmatp2_perp[temp[1][i][2]][temp[1][j][2]] + 1
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1) and (i == '1' or i == '7'):
										ggmatp2_b1[temp[1][i][2]][temp[1][j][2]]=ggmatp2_b1[temp[1][i][2]][temp[1][j][2]] + 1
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1) and (i == '2' or i == '8'):
										ggmatp2_b2[temp[1][i][2]][temp[1][j][2]]=ggmatp2_b2[temp[1][i][2]][temp[1][j][2]] + 1







	t1= Thread(target=matwrite,args=("ggmatp2_f.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2_f,'overwrite':overwrite})
	t2= Thread(target=matwrite,args=("ggmatp2_perp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2_perp,'overwrite':overwrite})
	t3= Thread(target=matwrite,args=("ggmatp2_b1.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2_b1,'overwrite':overwrite})
	t4= Thread(target=matwrite,args=("ggmatp2_b2.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2_b2,'overwrite':overwrite})
	
	
	t1.start()
	t2.start()
	t3.start()
	t4.start()
	
	
	t1.join()
	t2.join()
	t3.join()
	t4.join()

	

if __name__ == "__main__":
	import sys
	fileev5=sys.argv[1]
	overwrite=int(sys.argv[2])
	ev5toggmatmultqc(fileev5,dimy=5000,dimx=5000,overwrite=overwrite)
