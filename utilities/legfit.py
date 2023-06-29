#!/usr/bin/env python3
#Legendre Fit for Angular Distribution as a Function of Angle
#C. Wibisono 
#06/29/'23

#How to Use:
#./legfit [addata] A0_initial a2_initial a4_initial

import numpy as np
from scipy.optimize import curve_fit as cvt
from numpy.polynomial.legendre import legval as legpoly
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import sys
import os

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

def loaddata():
	global angle,intensity,gamma
	FILE=open(sys.argv[1])
	linefile=FILE.readlines()
	dim=len(linefile)
	angle=np.zeros(dim-1)
	intensity=np.zeros(dim-1)
	i=0
	for line in linefile:
		if line.find('#') == -1:
			liner=line.split()
			angle[i]=np.deg2rad(np.float32(liner[0]))
			intensity[i]=np.float32(liner[1])
			print(i,angle[i],intensity[i])
			i=i+1
		else:
			liner2=line.split()	
			gamma=liner2[0]
			print('Gamma Energy (keV):',gamma)
	FILE.close()


def legendrefunc(x,A0,a2,a4):
	P0=legpoly(np.cos(x),[1,0,0,0,0])
	P2=legpoly(np.cos(x),[0,0,1,0,0])
	P4=legpoly(np.cos(x),[0,0,0,0,1])
	return A0*(P0+a2*P2+a4*P4)

'''
def legendrefunc(x,A0,a2,a4):
	P0=1.
	P2=0.5*(3.*((np.cos(x))**2.)-1)
	P4=1./8.*(35.*((np.cos(x))**4.)-30.*((np.cos(x))**2.)+3.)
	return A0*(P0+a2*P2+a4*P4)
'''
def legendrefit():
	global popt,pcov
	p0=np.array([float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])])
	popt,pcov=cvt(legendrefunc,angle,intensity)
	print('A0:',popt[0],'a2:',popt[1],'a4:',popt[2])
	print('A0cov:',pcov[0],'a2cov:',pcov[1],'a4cov:',pcov[2])

def plot():
	x=np.arange(0,181,0.01)
	xtrf=np.deg2rad(x)
	fig,ax=plt.subplots()
	ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax.plot(np.cos(angle),intensity,'o',label='data')
	ax.plot(np.cos(xtrf),legendrefunc(xtrf,*popt),label='fit: A0=%5.3f, a2=%5.3f, a4=%5.3f' % tuple(popt))
	ax.legend()
	ax.set_xlabel(r'$\theta(rad)$',style='normal',fontweight='bold')
	ax.set_ylabel('Intensity(normalized)',style='normal',fontweight='bold')
	ax.set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	fig.suptitle("32P Angular Distribution\nClarion2-Trinity\nO16+O18 at 30 MeV")
	#ax.plot(np.cos(xtrf),legendrefunc(xtrf,11.09,-0.392,0.233),'m')
	plt.show()	

def main():
	loaddata()
	legendrefit()
	plot()

if __name__ == "__main__":
	main()
