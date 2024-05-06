#!/usr/bin/env python3

def linear(x,A,B):
	return A*x+B

def protoncalib(filecalib,ringnumber,calibresult):
	'''
	Proton Calibration Function:
	C. Wibisono
	Parameter(s):
	filecalib: calibration file consisting of trace integral and energy
	ringnumber: int 1 -5
	calibresult: calibration result file
	'''
	import matplotlib.pyplot as plt
	import numpy as np
	from scipy.optimize import curve_fit as cvt

	GAGG={}
	with open(filecalib,mode='r') as f:
		lines = f.readlines()
		for line in lines:
			if line.find('#') == -1:
				temp=list(line.split())
				if temp[0] in GAGG.keys():
					x=float(temp[1])
					y=float(temp[2])
					GAGG[temp[0]][0].append(x)
					GAGG[temp[0]][1].append(y)
				else:
					GAGG[temp[0]] =[[float(temp[1])],[float(temp[2])]]	
			else:
				pass
	popt=[]
	pcov=[]
	for i,j  in GAGG.items():
		popts,pcovs=cvt(linear, j[0], j[1])
		popt.append(popts)
		pcov.append(pcovs)
	
	for i in range(len(popt)):
		print(i+1,popt[i],pcov[i])

	xcoords=np.arange(0,4000,1)
	for i,j in GAGG.items():
		row=int(i)
		fig,ax=plt.subplots()
		ax.plot(GAGG[i][0],GAGG[i][1],'ro',label=i)
		ax.plot(xcoords,linear(xcoords,*popt[row-1]),color='b')
		ax.legend()
		plt.show()

	if ringnumber == 1:
		shift =0
	elif ringnumber == 2:
		shift = 8
	elif ringnumber == 3:
		shift = 18
	elif ringnumber == 4:
		shift = 32
	else:
		shift = 48
	with open(calibresult,mode='w') as fout:
		for i,j in GAGG.items():
			row=int(i)
			crystalnumber = row+shift
			fout.write(str(crystalnumber)+'\t'+str(popt[row-1][0])+'\t'+str(popt[row-1][1])+'\n')

if __name__ == "__main__":
	import sys
	fileproton=sys.argv[1]
	ringnumber = int(sys.argv[2])
	calibresult = sys.argv[3]
	protoncalib(fileproton,ringnumber,calibresult)		
