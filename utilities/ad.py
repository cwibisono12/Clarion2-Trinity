#!/usr/bin/env python3
#Legendre Fit for Angular Distribution as a Function of Angle
#Compare Theoretical Angular Distribution vs Mixing Ratio 
#C. Wibisono 
#06/29/'23
#07/11/'23 v2 with chisquare vs mixing ratio

#How to Use:
#./ad	[addata] A0_initial a2_initial a4_initial Ji Jf sigma

import numpy as np
from scipy.optimize import curve_fit as cvt
from numpy.polynomial.legendre import legval as legpoly
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import sys
import os
import adlib as ad

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

def loaddata():
	global angle,intensity,gamma,dim
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
	global popt,pcov,perr
	p0=np.array([float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4])])
	popt,pcov=cvt(legendrefunc,angle,intensity,p0)
	perr=np.sqrt(np.diag(pcov))
	print('A0:',popt[0],'a2:',popt[1],'a4:',popt[2])
	print('A0sderr:',perr[0],'a2sderr:',perr[1],'a4sderr:',perr[2])

def plot():
	x=np.arange(0,181,0.01)
	xtrf=np.deg2rad(x)
	fig,ax=plt.subplots()
	ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax.plot(np.cos(angle),intensity,'o',label='data')
	ax.plot(np.cos(xtrf),legendrefunc(xtrf,*popt),label='fit: A0=%5.3f, a2=%5.3f, a4=%5.3f' % tuple(popt))
	#ax.plot(np.cos(xtrf),legendrefunc(xtrf,11.09,-0.392,0.233),'m',label='AD')
	ax.legend()
	ax.set_xlabel(r'$cos(\theta)$',style='normal',fontweight='bold')
	ax.set_ylabel('Intensity(normalized)',style='normal',fontweight='bold')
	ax.set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	fig.suptitle("32P Angular Distribution\nClarion2-Trinity\nO16+O18 at 30 MeV")
	#plt.show()	


#Theoretical Intensity:
def theo(x,delta):
	global Ji,Jf,sigma
	Ji=int(sys.argv[5])
	Jf=int(sys.argv[6])
	sigma=float(sys.argv[7])
	B0=ad.Bk(Ji,0,sigma)
	B2=ad.Bk(Ji,2,sigma)
	B4=ad.Bk(Ji,4,sigma)
	R00=ad.Rk(0,Ji-Jf,Ji-Jf,Ji,Jf)
	R01=ad.Rk(0,Ji-Jf,Ji-Jf+1,Ji,Jf)
	R02=ad.Rk(0,Ji-Jf+1,Ji-Jf+1,Ji,Jf)
	R20=ad.Rk(2,Ji-Jf,Ji-Jf,Ji,Jf)
	R21=ad.Rk(2,Ji-Jf,Ji-Jf+1,Ji,Jf)
	R22=ad.Rk(2,Ji-Jf+1,Ji-Jf+1,Ji,Jf)
	R40=ad.Rk(4,Ji-Jf,Ji-Jf,Ji,Jf)
	R41=ad.Rk(4,Ji-Jf,Ji-Jf+1,Ji,Jf)
	R42=ad.Rk(4,Ji-Jf+1,Ji-Jf+1,Ji,Jf)
	Y0=B0*legpoly(np.cos(x),[1,0,0,0,0])*(R00+2.*R01*delta+R02*(delta**2.0))/(1.+(delta**2.))
	Y2=B2*legpoly(np.cos(x),[0,0,1,0,0])*(R20+2.*R21*delta+R22*(delta**2.0))/(1.+(delta**2.))
	Y4=B4*legpoly(np.cos(x),[0,0,0,0,1])*(R40+2.*R41*delta+R42*(delta**2.0))/(1.+(delta**2.))
	return Y0+Y2+Y4


def chisq(delta):
	chi=0.
	for l in range(0,dim-1,1):
		chi=chi+np.power(theo(angle[l],delta)-intensity[l],2.0)
	return chi

def plotchi():
	delta=np.arange(-3,3,0.001)
	mixratio=np.rad2deg(np.arctan(delta))
	fig2,ax2=plt.subplots()
	ax2.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax2.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax2.plot(mixratio,chisq(delta),label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf))
	ax2.legend()
	ax2.set_xlabel(r'$arctan(\delta)$',style='normal',fontweight='bold')
	ax2.set_ylabel('chisq',style='normal',fontweight='bold')
	ax2.set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	fig2.suptitle("32P Angular Distribution\nClarion2-Trinity\nO16+O18 at 30 MeV")
	#plt.show()	

def main():
	loaddata()
	legendrefit()
	plot()
	plotchi()
	plt.show()

if __name__ == "__main__":
	main()
