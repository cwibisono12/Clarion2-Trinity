#!/usr/bin/env python3
#Polarization Predict over Various Mixing Ratios
#C. Wibisono 08/08 '23


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import pollib2 as p

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

def pol(theta,Ji,Jf,mult,energy,Pexp,Perr,sigma):
	"""
	Interface to call polarization library pollib
	"""
	#P=np.zeros(181)
	mixratio=np.arange(-89,89,1.)
	delta=np.tan(np.deg2rad(mixratio))
	i=len(mixratio)
	P=np.zeros(i)
	for j in range(0,i,1):
		P[j]=p.polpredict(theta,delta[j],Ji,Jf,sigma,mult)
	'''
	fig,ax=plt.subplots()
	ax.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	ax.plot(mixratio,P,color='m',label='theory')
	ax.set_xlabel(r'$\tan^{-1}(\delta)$'+' '+'(degree)',style='normal',fontweight='bold')
	ax.set_ylabel('Polarization',style='normal',fontweight='bold')
	ax.axhline(y=Pexp,color='k',label='experiment')
	ax.axhspan(Pexp-Perr,Pexp+Perr,alpha=0.25,color='g')
	ax.set_ylim(-1,1)
	ax.legend(loc='upper right')
	fig.suptitle('Predicted Polarization\n'+'a2= '+str(round(a2,3))+' '+'a4= '+str(round(a4,3))+' '+r'$\theta$'+'='+' '+str(round(theta,1))+'\n'+'Ji:'+' '+str(Ji)+'---->'+'Jf:'+' '+str(Jf)+'\n'+'$E_\gamma$'+'='+' '+str(energy)+' '+'keV')
	plt.show()
	'''
	return P

if __name__ == "__main__":
	import sys
	theta=float(sys.argv[1])
	Ji=float(sys.argv[2])
	Jf=float(sys.argv[3])
	mult=int(sys.argv[4])
	energy=int(sys.argv[5])
	Pexp=float(sys.argv[6])
	Perr=float(sys.argv[7])
	sigma=float(sys.argv[8])
	pol(theta,Ji,Jf,mult,energy,Pexp,Perr,sigma)
