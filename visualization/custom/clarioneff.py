#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from scipy.optimize import curve_fit as cvt

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
plt.rcParams['figure.dpi'] = 200

def efffunc(energy,A,B,C):
	'''Relative Efficiency Function 
	for Clover Detectors'''
	return np.exp(A+B*np.log(energy)+C*np.power(np.log(energy),2.))

def getcoeff(fdata, fout):
	'''Function to Generate Coefficients of Relative Efficiency 
	for each Ring for Clover Detector
	Parameter(s):
	fdata: file input pointer object consisting of energy and counts for each Ring.
		see param/effdata_May2024.txt as an example
	fout: file output for printing parameter results for each Ring.
	'''
	Ring1counts = [] #48.24
	Ring2counts = [] #90
	Ring3counts = [] #131.75
	Ring4counts = [] #150
	xdata = []

	with open(fdata, mode='r') as f:
		lines = f.readlines()
		for line in lines:
			if line.find('#') == -1:
				temp = list(line.split())
				energy = float(temp[0])
				c1 = float(temp[1])
				c2 = float(temp[2])
				c3 = float(temp[3])
				c4 = float(temp[4])
				xdata.append(energy)
				Ring1counts.append(c1)
				Ring2counts.append(c2)
				Ring3counts.append(c3)
				Ring4counts.append(c4)

	poptRing1, pcovRing1 = cvt(efffunc, xdata, Ring1counts)
	poptRing2, pcovRing2 = cvt(efffunc, xdata, Ring2counts)
	poptRing3, pcovRing3 = cvt(efffunc, xdata, Ring3counts)
	poptRing4, pcovRing4 = cvt(efffunc, xdata, Ring4counts)
	
	with open(fout, mode = 'w') as f2:
		f2.write('#'+'Ring Angle'+'\t'+'A'+'\t'+'B'+'\t'+'C'+'\n')
		f2.write('48.24'+'\t'+str(poptRing1[0])+'\t'+str(poptRing1[1])+'\t'+str(poptRing1[2])+'\n')
		f2.write('90'+'\t'+str(poptRing2[0])+'\t'+str(poptRing2[1])+'\t'+str(poptRing2[2])+'\n')
		f2.write('131.75'+'\t'+str(poptRing3[0])+'\t'+str(poptRing3[1])+'\t'+str(poptRing3[2])+'\n')
		f2.write('150'+'\t'+str(poptRing4[0])+'\t'+str(poptRing4[1])+'\t'+str(poptRing4[2])+'\n')

	print('Fit Result: ')
	print('Ring1:',*poptRing1)
	print('Ring2:',*poptRing2)
	print('Ring3:',*poptRing3)
	print('Ring4:',*poptRing4)

	fig,ax=plt.subplots()
	ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

	xcoords=np.arange(100,2500,1)
	ax.plot(xcoords,efffunc(xcoords,*poptRing1),linewidth=0.85,color='b',label='48.24$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f' % tuple(poptRing1))
	ax.plot(xcoords,efffunc(xcoords,*poptRing2),linewidth=0.85,color='k',label='90$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f' % tuple(poptRing2))
	ax.plot(xcoords,efffunc(xcoords,*poptRing3),linewidth=0.85,color='g',label='131.75$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f' % tuple(poptRing3))
	ax.plot(xcoords,efffunc(xcoords,*poptRing4),linewidth=0.85,color='r',label='150$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f' % tuple(poptRing4))

	ax.plot(xdata,Ring1counts,'bo')
	ax.plot(xdata,Ring2counts,'ko')
	ax.plot(xdata,Ring3counts,'go')
	ax.plot(xdata,Ring4counts,'ro')
	ax.legend()
	ax.set_xlabel(r'E$_{\gamma}$ (keV)',fontweight='bold')
	ax.set_ylabel(r'$\epsilon_{eff}$ (arb)',fontweight='bold')
	ax.text(1500,2.0,'$\epsilon_{eff} = e^{A+B\log(E)+C(\log(E))^{2}}$')
	ax.set_xlim(100,2500)
	plt.show()

def clarioneff(fparam, energy):
	'''Function to deduce the relative efficiency given the energy for each Ring of Clover array
	Parameter(s):
	fparam: parameter file pointer object
	energy: gamma energy
	Return(s):
	effRing: (dict) Relative efficiency at a given energy for each Ring
	'''
	coeffRing = {}
	effRing = {}
	with open(fparam, mode='r') as f:
		lines = f.readlines()
		for line in lines:
			if line.find('#') == -1:
				temp = list(line.split())
				coeffRing[temp[0]] = [float(temp[1]),float(temp[2]),float(temp[3])] 		

	for i, j in coeffRing.items():
		effRing[i] = efffunc(energy,j[0],j[1],j[2])

	return effRing


if __name__ == "__main__":
	import sys
	filein = sys.argv[1]
	fileout = sys.argv[2]
	getcoeff(filein,fileout)
