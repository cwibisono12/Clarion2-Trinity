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
	
def efffunc2(energy,A,B,C,D):
	'''Another variant of Relative Efficiency Function 
	for Clover Detectors'''
	return np.exp(A+B*np.log(energy)+C*np.power(np.log(energy),2.)+D/(energy**2.))

def getcoeff_total(fdata, fout, norma, normb):
	'''Function to Generate Coefficients of Relative Efficiency 
	for the entire Arrays of HPGe Clover detectors.
	Parameter(s):
	fdata: file input pointer object consisting of energy and counts.
		see param/effdata_Aug2024b.txt as an example
	fout: file output for printing parameter results.
	norma: normalization factor applied to the first source.
	normb: normalization factor applied to the second souce.
	'''
	
	data = {}

	with open(fdata, mode='r') as f:
		lines = f.readlines()
		for line in lines:
			if line.find('#') == -1:
				temp = list(line.split())
				source = temp[0]
				if source in data.keys():
					data[source][0].append(float(temp[1])) #energy
					data[source][1].append(float(temp[2])) #relative eff
				else:
					data[source] = [[float(temp[1])],[float(temp[2])]]
				
	fig,ax=plt.subplots()
	ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	xdata = []
	ydata = []
	for name in data.keys():
		if name == '56Co':
			colorname = 'r+'
			dim = len(data[name][1])
			for k in range(dim):
				data[name][1][k] = data[name][1][k]*normb
				xdata.append(data[name][0][k])
				ydata.append(data[name][1][k])
		else:
			colorname = 'bo'
			dim = len(data[name][1])
			for k in range(dim):
				data[name][1][k] = data[name][1][k]*norma
				xdata.append(data[name][0][k])
				ydata.append(data[name][1][k])

		ax.plot(data[name][0],data[name][1],colorname,label=name)
	
	
	
	popt, pcov = cvt(efffunc2, xdata, ydata)
	
	with open(fout, mode = 'w') as f2:
		f2.write('#'+'Relative Eff'+'\t'+'A'+'\t'+'B'+'\t'+'C'+'\t'+'D'+'\n')
		f2.write('all'+'\t'+str(popt[0])+'\t'+str(popt[1])+'\t'+str(popt[2])+'\t'+str(popt[3])+'\n')

	print('Fit Result: ')
	print('popt:',*popt)


	xcoords=np.arange(1,4500,1)
	ax.plot(xcoords,efffunc2(xcoords,*popt),linewidth=0.85,color='b',label='fit: A=%5.3f, B=%5.3f, C=%5.3f, D=%5.3f' % tuple(popt))
	
	ax.legend()
	ax.set_xlabel(r'E$_{\gamma}$ (keV)',fontweight='bold')
	ax.set_ylabel(r'$\epsilon_{eff}$ (arb)',fontweight='bold')
	#ax.text(1500,1.0,'$\epsilon_{eff} = 10^{A}(E^{B})(E^{2})^{C}(10^{D/E^{2}})$')
	ax.text(1500,1.0,'$\epsilon_{eff} = e^{A+B\log(E)+C(\log(E))^{2}+D/E^{2}}$')
	ax.set_xlim(1,4500)
	plt.show()


def getcoeff(fdata, fout, norma, normb):
	'''Function to Generate Coefficients of Relative Efficiency 
	for each Ring for Clover Detector
	Parameter(s):
	fdata: file input pointer object consisting of energy and counts for each Ring.
		see param/effdata_Aug2024_ring.txt as an example
	fout: file output for printing parameter results for each Ring.
	'''
	Ring1counts = [] #48.24
	Ring2counts = [] #90
	Ring3counts = [] #131.75
	Ring4counts = [] #150
	xdata = []

	data = {}

	with open(fdata, mode='r') as f:
		lines = f.readlines()
		for line in lines:
			if line.find('#') == -1:
				temp = list(line.split())
				source = temp[0]
				if source in data.keys():
					data[source][0].append(float(temp[1])) #energy
					data[source][1].append(float(temp[2])) #intensity at 48.24
					data[source][2].append(float(temp[3])) #intensity at 90
					data[source][3].append(float(temp[4])) #intensity at 131.75
					data[source][4].append(float(temp[5])) #intensity at 150
				else:
					data[source] = [[float(temp[1])],[float(temp[2])],[float(temp[3])],[float(temp[4])],[float(temp[5])]]

	fig,ax=plt.subplots()
	ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	ax.yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	
	for name in data.keys():
		if name == '56Co':
			colorname = ['r+','r*','ro','r>']
			dim = len(data[name][1])
			for k in range(dim):
				data[name][1][k] = data[name][1][k]*normb
				data[name][2][k] = data[name][2][k]*normb
				data[name][3][k] = data[name][3][k]*normb
				data[name][4][k] = data[name][4][k]*normb
				xdata.append(data[name][0][k])
				Ring1counts.append(data[name][1][k])
				Ring2counts.append(data[name][2][k])
				Ring3counts.append(data[name][3][k])
				Ring4counts.append(data[name][4][k])
		else:
			colorname = ['b+','b*','bo','b>']
			dim = len(data[name][1])
			for k in range(dim):
				data[name][1][k] = data[name][1][k]*norma
				data[name][2][k] = data[name][2][k]*norma
				data[name][3][k] = data[name][3][k]*norma
				data[name][4][k] = data[name][4][k]*norma
				xdata.append(data[name][0][k])
				Ring1counts.append(data[name][1][k])
				Ring2counts.append(data[name][2][k])
				Ring3counts.append(data[name][3][k])
				Ring4counts.append(data[name][4][k])

		#ax.plot(data[name][0],data[name][1],colorname[0],label=name+' '+'48.24')
		#ax.plot(data[name][0],data[name][2],colorname[1],label=name+' '+'90')
		#ax.plot(data[name][0],data[name][3],colorname[2],label=name+' '+'131.75')
		#ax.plot(data[name][0],data[name][4],colorname[3],label=name+' '+'150')
		
	
	poptRing1, pcovRing1 = cvt(efffunc2, xdata, Ring1counts)
	poptRing2, pcovRing2 = cvt(efffunc2, xdata, Ring2counts)
	poptRing3, pcovRing3 = cvt(efffunc2, xdata, Ring3counts)
	poptRing4, pcovRing4 = cvt(efffunc2, xdata, Ring4counts)
	
	with open(fout, mode = 'w') as f2:
		f2.write('#'+'Ring Angle'+'\t'+'A'+'\t'+'B'+'\t'+'C'+'\t'+'D'+'\n')
		f2.write('48.24'+'\t'+str(poptRing1[0])+'\t'+str(poptRing1[1])+'\t'+str(poptRing1[2])+'\t'+str(poptRing1[3])+'\n')
		f2.write('90'+'\t'+str(poptRing2[0])+'\t'+str(poptRing2[1])+'\t'+str(poptRing2[2])+'\t'+str(poptRing2[3])+'\n')
		f2.write('131.75'+'\t'+str(poptRing3[0])+'\t'+str(poptRing3[1])+'\t'+str(poptRing3[2])+'\t'+str(poptRing3[3])+'\n')
		f2.write('150'+'\t'+str(poptRing4[0])+'\t'+str(poptRing4[1])+'\t'+str(poptRing4[2])+'\t'+str(poptRing4[3])+'\n')

	print('Fit Result: ')
	print('Ring1:',*poptRing1)
	print('Ring2:',*poptRing2)
	print('Ring3:',*poptRing3)
	print('Ring4:',*poptRing4)


	xcoords=np.arange(1,4500,1)
	ax.plot(xcoords,efffunc2(xcoords,*poptRing1),linewidth=0.85,color='b',label='48.24$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f, D=%5.3f' % tuple(poptRing1))
	ax.plot(xcoords,efffunc2(xcoords,*poptRing2),linewidth=0.85,color='k',label='90$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f, D=%5.3f' % tuple(poptRing2))
	ax.plot(xcoords,efffunc2(xcoords,*poptRing3),linewidth=0.85,color='g',label='131.75$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f, D=%5.3f' % tuple(poptRing3))
	ax.plot(xcoords,efffunc2(xcoords,*poptRing4),linewidth=0.85,color='r',label='150$^{0}$, fit: A=%5.3f, B=%5.3f, C=%5.3f, D=%5.3f' % tuple(poptRing4))

	ax.plot(xdata,Ring1counts,'bo')
	ax.plot(xdata,Ring2counts,'ko')
	ax.plot(xdata,Ring3counts,'go')
	ax.plot(xdata,Ring4counts,'ro')
	ax.legend()
	ax.set_xlabel(r'E$_{\gamma}$ (keV)',fontweight='bold')
	ax.set_ylabel(r'$\epsilon_{eff}$ (arb)',fontweight='bold')
	ax.text(1500,4.0,'$\epsilon_{eff} = e^{A+B\log(E)+C(\log(E))^{2}+D/E^{2}}$')
	ax.set_xlim(1,4500)
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
				coeffRing[temp[0]] = [float(temp[1]),float(temp[2]),float(temp[3]),float(temp[4])] 		

	for i, j in coeffRing.items():
		effRing[i] = efffunc2(energy,j[0],j[1],j[2],j[3])

	return effRing


if __name__ == "__main__":
	import sys
	filein = sys.argv[1]
	fileout = sys.argv[2]
	#norma = float(sys.argv[2])
	#normb = float(sys.argv[3])
	#eff = getcoeff_total(filein,fileout,  1.0, 0.34)
	getcoeff(filein,fileout,1.0,0.345)
	#print(eff)
