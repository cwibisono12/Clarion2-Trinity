#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread

def ev5toggmatmultq(fpr,*,dimy=5000,dimx=5000,overwrite=0):
	'''
	Function to Generate ggmatrix from ev5 file
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension
	dimx: x dimension
	overwrite: overwrite flag 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	ggmatp2=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp2_nd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatp2_dd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#Non-Prompt:
	ggmatpa=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatpa_nd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmatpa_dd=np.ndarray(shape=(dimy,dimx),dtype=np.int32)

	#tdiff:
	ggmat_tdiffp2=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	ggmat_tdiffpa=np.ndarray(shape=(dimy,dimx),dtype=np.int32)

	for i in range(dimy):
		for j in range(dimx):
			ggmatp2[i][j]=0
			ggmatp2_nd[i][j]=0
			ggmatp2_dd[i][j]=0
			ggmatpa[i][j]=0
			ggmatpa_nd[i][j]=0
			ggmatpa_dd[i][j]=0
			ggmat_tdiffp2[i][j]=0
			ggmat_tdiffpa[i][j]=0

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
								if (tdiff > 0 and tdiff <5000) and (temp[1][i][2] >=0 and temp[1][i][2] <5000)	and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
									if temp[1][i][2] < temp[1][j][2]:
										ggmat_tdiffp2[temp[1][i][2]][tdiff] = ggmat_tdiffp2[temp[1][i][2]][tdiff] + 1	
									else:											
										ggmat_tdiffp2[temp[1][j][2]][tdiff] = ggmat_tdiffp2[temp[1][i][2]][tdiff] + 1	



								if tdiff >= 1988 and tdiff <=2012:
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1):
										ggmatp2[temp[1][i][2]][temp[1][j][2]]=ggmatp2[temp[1][i][2]][temp[1][j][2]] + 1	
									if (temp[1][i][3] >= 0 and temp[1][i][3] <5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1):
										ggmatp2_nd[temp[1][i][3]][temp[1][j][2]]=ggmatp2_nd[temp[1][i][3]][temp[1][j][2]] + 1
									if (temp[1][i][3] >= 0 and temp[1][i][3] < 5000) and (temp[1][j][3] >=0 and temp[1][j][3] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1):
										ggmatp2_dd[temp[1][i][3]][temp[1][j][3]]=ggmatp2_dd[temp[1][i][3]][temp[1][j][3]] + 1


				if temp[0][0] == 1 and temp[0][1] == 1: #1p-1a gated
					gmult=temp[0][2]
					for i in temp[1].keys():
						for j in temp[1].keys():
							if i!=j :
								tdiff=temp[1][i][6] - temp[1][j][6] + 2000
								if (tdiff > 0 and tdiff <5000) and (temp[1][i][2] >=0 and temp[1][i][2] <5000)	and (temp[1][j][2] >=0 and temp[1][j][2] <5000):
									if temp[1][i][2] < temp[1][j][2]:
										ggmat_tdiffpa[temp[1][i][2]][tdiff] = ggmat_tdiffpa[temp[1][i][2]][tdiff] + 1	
									else:											
										ggmat_tdiffpa[temp[1][j][2]][tdiff] = ggmat_tdiffpa[temp[1][i][2]][tdiff] + 1	



								if tdiff >= 1988 and tdiff <=2012:
									if (temp[1][i][2] >= 0 and temp[1][i][2] < 5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1):
										ggmatpa[temp[1][i][2]][temp[1][j][2]]=ggmatpa[temp[1][i][2]][temp[1][j][2]] + 1	
									if (temp[1][i][3] >= 0 and temp[1][i][3] <5000) and (temp[1][j][2] >=0 and temp[1][j][2] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1):
										ggmatpa_nd[temp[1][i][3]][temp[1][j][2]]=ggmatpa_nd[temp[1][i][3]][temp[1][j][2]] + 1
									if (temp[1][i][3] >= 0 and temp[1][i][3] < 5000) and (temp[1][j][3] >=0 and temp[1][j][3] <5000) and (temp[1][i][4] == 1 and temp[1][j][4] == 1):
										ggmatpa_dd[temp[1][i][3]][temp[1][j][3]]=ggmatpa_dd[temp[1][i][3]][temp[1][j][3]] + 1





	t1= Thread(target=matwrite,args=("ggmat_tdiffp2.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmat_tdiffp2,'overwrite':overwrite})
	t2= Thread(target=matwrite,args=("ggmatp2.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2,'overwrite':overwrite})
	t3= Thread(target=matwrite,args=("ggmatp2_nd.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2_nd,'overwrite':overwrite})
	t4= Thread(target=matwrite,args=("ggmatp2_dd.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatp2_dd,'overwrite':overwrite})
	
	t1b= Thread(target=matwrite,args=("ggmat_tdiffpa.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmat_tdiffpa,'overwrite':overwrite})
	t2b= Thread(target=matwrite,args=("ggmatpa.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatpa,'overwrite':overwrite})
	t3b= Thread(target=matwrite,args=("ggmatpa_nd.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatpa_nd,'overwrite':overwrite})
	t4b= Thread(target=matwrite,args=("ggmatpa_dd.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':ggmatpa_dd,'overwrite':overwrite})
	
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
	fileev5=sys.argv[1]
	overwrite=int(sys.argv[2])
	ev5toggmatmultq(fileev5,dimy=5000,dimx=5000,overwrite=overwrite)
