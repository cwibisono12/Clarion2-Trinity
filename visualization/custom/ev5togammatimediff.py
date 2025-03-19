#!/usr/bin/env python3

import pxi16parser as p
from parsegen import matwrite
import numpy as np
from threading import Thread 

def ev5togammatimediff(fpr,*,dimy=5000,dimx=5000,overwrite=0):
	'''
	Function to Generate Gamma Energy vs TimeDifference (measured with respect to particle emission) from ev5 file
	Parameter(s):
	fpr: ev5 file pointer object
	dimy: y dimension
	dimx: x dimension
	overwrite: 1 (to overwrite) or 0 (to append)
	'''
	#Prompt:
	p2_gammatdp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_gammatdp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	p2_gammatdp_dopp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	pa_gammatdp_dopp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	#p2_Exgammap_dopp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	#pa_Exgammap_dopp = np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	
	
	for i in range(dimy):
		for j in range(dimx):
			p2_gammatdp[i][j]=0
			pa_gammatdp[i][j]=0
			
			
			p2_gammatdp_dopp[i][j]=0
			pa_gammatdp_dopp[i][j]=0
			
			

	with open(fpr,mode='rb') as f:
		while(1):
			temp=p.ev5read(f)
			if temp == -1:
				break
			else:
				if (temp[0][0] == 2 and temp[0][1] == 0): #2p-gated
					#temp[0][3] = temp[0][3] - 180#to include energy loss thorugh the target (1.8 MeV for 30 MeV and 1.58 for 32MeV)
					#temp[0][3] = temp[0][3] - 250#to include energy loss thorugh the target (2.8 MeV for 32 MeV O18+O18)
					
					for i in temp[1].keys():
						timediff = temp[1][i][6]
						energy = temp[1][i][2]
						edopp = temp[1][i][3]
						if ((energy >=0 and energy < 5000) and (timediff>= 0 and timediff < 5000) and (temp[1][i][4] == 1 and temp[1][i][5]==0)):
							p2_gammatdp[energy][timediff] = p2_gammatdp[energy][timediff] + 1
						if ((edopp >=0 and edopp < 5000) and (timediff>= 0 and timediff < 5000) and (temp[1][i][4] == 1 and temp[1][i][5]==0)):
							p2_gammatdp_dopp[edopp][timediff] = p2_gammatdp_dopp[edopp][timediff] + 1

				if (temp[0][0] == 1 and temp[0][1] == 1): #1p-1agated
					#temp[0][3] = temp[0][3] - 395 #to inlcude energy loss through the target 
					i#temp[0][3] = temp[0][3] - 445
					for i in temp[1].keys():
						timediff = temp[1][i][6]
						energy = temp[1][i][2]
						edopp = temp[1][i][3]
						if ((energy >=0 and energy <5000) and (timediff>=0 and timediff<5000) and (temp[1][i][4] == 1 and temp[1][i][5] ==0)):
							pa_gammatdp[energy][timediff] = pa_gammatdp[energy][timediff] + 1
						if ((edopp >=0 and edopp <5000) and (timediff>=0 and timediff<5000) and (temp[1][i][4] == 1 and temp[1][i][5] ==0)):
							pa_gammatdp_dopp[edopp][timediff] = pa_gammatdp_dopp[edopp][timediff] + 1
						

			
	t1 = Thread(target=matwrite,args=("p2_gammatdp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_gammatdp,'overwrite':overwrite})
	t2 = Thread(target=matwrite,args=("pa_gammatdp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_gammatdp,'overwrite':overwrite})
	t3 = Thread(target=matwrite,args=("p2_gammatdp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':p2_gammatdp_dopp,'overwrite':overwrite})
	t4 = Thread(target=matwrite,args=("pa_gammatdp_dopp.spn2",),kwargs={'dimy':dimy,'dimx':dimx,'arr':pa_gammatdp_dopp,'overwrite':overwrite})
	
	
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

	ev5togammatimediff(fileev5,dimy=5000,dimx=5000,overwrite=overwrite)
