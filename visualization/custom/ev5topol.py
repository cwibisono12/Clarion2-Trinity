#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread 

def ev5topol(fpr,*,dimy=2,dimx=8192,overwrite=0):
	'''
	Function to Generate Linear Polarization Spectra from ev5 file
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension 0->Paralel, 1->Perpendicular
	dimx: x dimension
	overwrite: 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	p2_polp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_polp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p2_polp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_polp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#Non-Prompt:
	p2_polnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_polnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	p2_polnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_polnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	for i in range(dimy):
		for j in range(dimx):
			p2_polp[i][j]=0
			pa_polp[i][j]=0
			
			p2_polp_dopp[i][j]=0
			pa_polp_dopp[i][j]=0
			
			p2_polnp[i][j]=0
			pa_polnp[i][j]=0
			
			p2_polnp_dopp[i][j]=0
			pa_polnp_dopp[i][j]=0

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if temp[0][0] == 2 and temp[0][1] == 0: #2p-gated
					for i in temp[1].keys():
						#Prompt:
						if (temp[1][i][4] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
							#Parallel Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
								paralel = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (paralel >= 0 and paralel < 8192):
									p2_polp[0][paralel] = p2_polp[0][paralel] + 1
								if (edopp >= 0 and edopp < 8192):
									p2_polp_dopp[0][edopp] = p2_polp_dopp[0][edopp] + 1

							#Perpendicular Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
								perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (perpendicular >= 0 and perpendicular < 8192):
									p2_polp[1][perpendicular] = p2_polp[1][perpendicular] + 1
								if (edopp >= 0 and edopp < 8192):
									p2_polp_dopp[1][edopp] = p2_polp_dopp[1][edopp] + 1
						#Non-Prompt:
						if (temp[1][i][5] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
							#Parallel Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
								paralel = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (paralel >= 0 and paralel < 8192):
									p2_polnp[0][paralel] = p2_polnp[0][paralel] + 1
								if (edopp >= 0 and edopp < 8192):
									p2_polnp_dopp[0][edopp] = p2_polnp_dopp[0][edopp] + 1

							#Perpendicular Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
								perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (perpendicular >= 0 and perpendicular < 8192):
									p2_polnp[1][perpendicular] = p2_polnp[1][perpendicular] + 1
								if (edopp >= 0 and edopp < 8192):
									p2_polnp_dopp[1][edopp] = p2_polnp_dopp[1][edopp] + 1

				if temp[0][0] == 1 and temp[0][1] == 1: #1p-1a gated
					for i in temp[1].keys():
						#Prompt:
						if (temp[1][i][4] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
							#Parallel Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
								paralel = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (paralel >= 0 and paralel < 8192):
									pa_polp[0][paralel] = pa_polp[0][paralel] + 1
								if (edopp >= 0 and edopp < 8192):
									pa_polp_dopp[0][edopp] = pa_polp_dopp[0][edopp] + 1

							#Perpendicular Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
								perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (perpendicular >= 0 and perpendicular < 8192):
									pa_polp[1][perpendicular] = pa_polp[1][perpendicular] + 1
								if (edopp >= 0 and edopp < 8192):
									pa_polp_dopp[1][edopp] = pa_polp_dopp[1][edopp] + 1
						#Non-Prompt:
						if (temp[1][i][5] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
							#Parallel Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
								paralel = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (paralel >= 0 and paralel < 8192):
									pa_polnp[0][paralel] = pa_polnp[0][paralel] + 1
								if (edopp >= 0 and edopp < 8192):
									pa_polnp_dopp[0][edopp] = pa_polnp_dopp[0][edopp] + 1

							#Perpendicular Spectra:
							if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
								perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
								edopp = temp[1][i][3] 
								if (perpendicular >= 0 and perpendicular < 8192):
									pa_polnp[1][perpendicular] = pa_polnp[1][perpendicular] + 1
								if (edopp >= 0 and edopp < 8192):
									pa_polnp_dopp[1][edopp] = pa_polnp_dopp[1][edopp] + 1


			
	t1 = Thread(target=matwrite,args=("p2_polp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_polp,'overwrite':overwrite})
	t2 = Thread(target=matwrite,args=("p2_polp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_polp_dopp,'overwrite':overwrite})
	t3 = Thread(target=matwrite,args=("pa_polp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_polp,'overwrite':overwrite})
	t4 = Thread(target=matwrite,args=("pa_polp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_polp_dopp,'overwrite':overwrite})
	
	t1b = Thread(target=matwrite,args=("p2_polnp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_polnp,'overwrite':overwrite})
	t2b = Thread(target=matwrite,args=("p2_polnp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_polnp_dopp,'overwrite':overwrite})
	t3b = Thread(target=matwrite,args=("pa_polnp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_polnp,'overwrite':overwrite})
	t4b = Thread(target=matwrite,args=("pa_polnp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_polnp_dopp,'overwrite':overwrite})
	
	t1.start()
	t2.start()
	t3.start()
	t4.start()

	t1b.start()
	t2b.start()
	t3b.start()
	t4b.start()

	t1.join()
	t2.join()
	t3.join()
	t4.join()

	t1b.join()
	t2b.join()
	t3b.join()
	t4b.join()


if __name__ == "__main__":
	import sys
	fileev5 = sys.argv[1]
	overwrite = int(sys.argv[2])

	ev5topol(fileev5,dimy=16,dimx=8192,overwrite=overwrite)
