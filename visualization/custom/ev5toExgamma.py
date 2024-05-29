#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread 

def ev5toExgamma(fpr,*,dimy=2500,dimx=8192,overwrite=0):
	'''
	Function to Generate Recoil Excitation Energy-Gamma Energies from ev5 file
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension
	dimx: x dimension
	overwrite: 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	p2_Exgammap = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_Exgammap = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p2_Exgammap_dopp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_Exgammap_dopp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	
	for i in range(dimy):
		for j in range(dimx):
			p2_Exgammap[i][j]=0
			pa_Exgammap[i][j]=0
			
			
			p2_Exgammap_dopp[i][j]=0
			pa_Exgammap_dopp[i][j]=0
			
			

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if (temp[0][0] == 2 and temp[0][1] == 0) and (temp[0][3] > 0 and temp[0][3] <= 2500): #2p-gated
					for i in temp[1].keys():
						if (temp[1][i][2] >=0 and temp[1][i][2] < 8192) and (temp[1][i][4] == 1):
							p2_Exgammap[temp[0][3]][temp[1][i][2]] = p2_Exgammap[temp[0][3]][temp[1][i][2]] + 1
						if (temp[1][i][3] >=0 and temp[1][i][3] < 8192) and (temp[1][i][4] == 1):
							p2_Exgammap_dopp[temp[0][3]][temp[1][i][3]]=p2_Exgammap_dopp[temp[0][3]][temp[1][i][3]] + 1

				if (temp[0][0] == 1 and temp[0][1] == 1) and (temp[0][3] > 0 and temp[0][3] <= 2500): #1p-1agated
					for i in temp[1].keys():
						if (temp[1][i][2] >=0 and temp[1][i][2] <8192) and (temp[1][i][4] == 1):
							pa_Exgammap[temp[0][3]][temp[1][i][2]] = pa_Exgammap[temp[0][3]][temp[1][i][2]]+1
						
						if (temp[1][i][3] >=0 and temp[1][i][3] <8192) and (temp[1][i][4] == 1):
							pa_Exgammap_dopp[temp[0][3]][temp[1][i][3]]=pa_Exgammap_dopp[temp[0][3]][temp[1][i][3]]+1

			
	t1 = Thread(target=matwrite,args=("p2_Exgammap.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_Exgammap,'overwrite':overwrite})
	t2 = Thread(target=matwrite,args=("p2_Exgammap_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_Exgammap_dopp,'overwrite':overwrite})
	t3 = Thread(target=matwrite,args=("pa_Exgammap.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_Exgammap,'overwrite':overwrite})
	t4 = Thread(target=matwrite,args=("pa_Exgammap_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_Exgammap_dopp,'overwrite':overwrite})
	
	
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
	fileev5 = sys.argv[1]
	overwrite = int(sys.argv[2])

	ev5toExgamma(fileev5,dimy=2500,dimx=8192,overwrite=overwrite)
