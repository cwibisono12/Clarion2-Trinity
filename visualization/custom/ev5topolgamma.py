#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread 

def ev5topolg(fpr,xlow,xup,xlowbkg,xupbkg,*,dimy=2,dimx=8192,overwrite=0):
	'''
	Function to Generate Linear Polarization Spectra from ev5 file gated on particle and gamma
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension 0->Paralel, 1->Perpendicular
	dimx: x dimension
	overwrite: 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	p1_polp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polp_dopp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polp_dopp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#Non-Prompt:
	p1_polnp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polnp_dopp=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polnp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polnp_dopp_p=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polnp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p1_polnp_dopp_np=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	for i in range(dimy):
		for j in range(dimx):
			p1_polp[i][j]=0
			p1_polp_dopp[i][j]=0
			p1_polp_p[i][j]=0
			p1_polp_dopp_p[i][j]=0
			p1_polp_np[i][j]=0
			p1_polp_dopp_np[i][j]=0
			p1_polnp[i][j]=0
			p1_polnp_dopp[i][j]=0
			p1_polnp_p[i][j]=0
			p1_polnp_dopp_p[i][j]=0
			p1_polnp_np[i][j]=0
			p1_polnp_dopp_np[i][j]=0

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if temp[0][0] == 1 and temp[0][1] == 0: #1p-gated
					for j in temp[1].keys():
						'''Prompt Region'''
						if((temp[1][j][2] >=xlow and temp[1][j][2] <= xup) and temp[1][j][4] == 1): 

							for i in temp[1].keys():
								if (i!=j):
									#Prompt:
									if (temp[1][i][4] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
										#Parallel Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
											paralel = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (paralel >= 0 and paralel < 8192):
												p1_polp_p[0][paralel] = p1_polp_p[0][paralel] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polp_dopp_p[0][edopp] = p1_polp_dopp_p[0][edopp] + 1

										#Perpendicular Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
											perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (perpendicular >= 0 and perpendicular < 8192):
												p1_polp_p[1][perpendicular] = p1_polp_p[1][perpendicular] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polp_dopp_p[1][edopp] = p1_polp_dopp_p[1][edopp] + 1
						
									#Non-Prompt:
									if (temp[1][i][5] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
									#Parallel Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
											paralel = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (paralel >= 0 and paralel < 8192):
												p1_polnp_p[0][paralel] = p1_polnp_p[0][paralel] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polnp_dopp_p[0][edopp] = p1_polnp_dopp_p[0][edopp] + 1

									#Perpendicular Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
											perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (perpendicular >= 0 and perpendicular < 8192):
												p1_polnp_p[1][perpendicular] = p1_polnp_p[1][perpendicular] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polnp_dopp_p[1][edopp] = p1_polnp_dopp_p[1][edopp] + 1
						
						'''Background Region'''
						if((temp[1][j][2] >=xlowbkg and temp[1][j][2] <= xupbkg) and temp[1][j][4] == 1): 

							for i in temp[1].keys():
								if (i!=j):
									#Prompt:
									if (temp[1][i][4] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
										#Parallel Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
											paralel = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (paralel >= 0 and paralel < 8192):
												p1_polp_np[0][paralel] = p1_polp_np[0][paralel] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polp_dopp_np[0][edopp] = p1_polp_dopp_np[0][edopp] + 1

										#Perpendicular Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
											perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (perpendicular >= 0 and perpendicular < 8192):
												p1_polp_np[1][perpendicular] = p1_polp_np[1][perpendicular] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polp_dopp_np[1][edopp] = p1_polp_dopp_np[1][edopp] + 1
						
									#Non-Prompt:
									if (temp[1][i][5] == 1 and temp[1][i][7] == 2) and (i == '3' or i == '4' or i== '9' or i == '10'):
									#Parallel Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 3) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 2):
											paralel = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (paralel >= 0 and paralel < 8192):
												p1_polnp_np[0][paralel] = p1_polnp_np[0][paralel] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polnp_dopp_np[0][edopp] = p1_polnp_dopp_np[0][edopp] + 1

									#Perpendicular Spectra:
										if (temp[1][i][8][0] == 1 and temp[1][i][8][1] == 2) or (temp[1][i][8][0] == 2 and temp[1][i][8][1] == 1) or (temp[1][i][8][0] == 3 and temp[1][i][8][1] == 4) or (temp[1][i][8][0] == 4 and temp[1][i][8][1] == 3):
											perpendicular = temp[1][i][9][0] + temp[1][i][9][1]
											edopp = temp[1][i][3] 
											if (perpendicular >= 0 and perpendicular < 8192):
												p1_polnp_np[1][perpendicular] = p1_polnp_np[1][perpendicular] + 1
											if (edopp >= 0 and edopp < 8192):
												p1_polnp_dopp_np[1][edopp] = p1_polnp_dopp_np[1][edopp] + 1



	for k1 in range(dimy):
		for k2 in range(dimx):
			p1_polp[k1][k2] = p1_polp_p[k1][k2] - p1_polnp_p[k1][k2]
			p1_polp_dopp[k1][k2] = p1_polp_dopp_p[k1][k2] - p1_polnp_dopp_p[k1][k2]
			p1_polnp[k1][k2] = p1_polp_np[k1][k2] - p1_polnp_np[k1][k2]
			p1_polnp_dopp[k1][k2] = p1_polp_dopp_np[k1][k2] - p1_polnp_dopp_np[k1][k2]
			
	t1 = Thread(target=matwrite,args=("p1_polp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_polp,'overwrite':overwrite})
	t2 = Thread(target=matwrite,args=("p1_polp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_polp_dopp,'overwrite':overwrite})
	
	t1b = Thread(target=matwrite,args=("p1_polnp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_polnp,'overwrite':overwrite})
	t2b = Thread(target=matwrite,args=("p1_polnp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p1_polnp_dopp,'overwrite':overwrite})
	
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

	ev5topolg(fileev5,xlow,xup,xlowbkg,xupbkg,dimy=2,dimx=8192,overwrite=overwrite)
