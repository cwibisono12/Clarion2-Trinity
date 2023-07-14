#!/usr/bin/env python3
#Legendre Fit for Angular Distribution as a Function of Angle
#Compare Theoretical Angular Distribution vs Mixing Ratio 
#C. Wibisono 
#06/29/'23
#07/11/'23 v2 with chisquare vs mixing ratio

#How to Use:
#./ad [addata] [list of Ji/list of m] A0_initial a2_initial a4_initial Jf [option] Ji/m
#option == 1 for vary Ji fix m
#option == 2 for vary m fix Ji

import numpy as np
from scipy.optimize import curve_fit as cvt
from numpy.polynomial.legendre import legval as legpoly
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import sys
import os
import adlib as ad
import datetime as dt

plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']


option = int(sys.argv[7])

def loaddata():
	global angle,intensity,error,gamma,dim
	FILE=open(sys.argv[1])
	linefile=FILE.readlines()
	dim=len(linefile)
	angle=np.zeros(dim-1)
	intensity=np.zeros(dim-1)
	error=np.zeros(dim-1)
	i=0
	for line in linefile:
		if line.find('#') == -1:
			liner=line.split()
			angle[i]=np.deg2rad(np.float32(liner[0]))
			intensity[i]=np.float32(liner[1])
			error[i]=np.float32(liner[2])
			print(i,angle[i],intensity[i],error[i])
			i=i+1
		else:
			liner2=line.split()	
			gamma=liner2[0]
			print('Gamma Energy (keV):',gamma)
	FILE.close()

def loadJi():
	global Jinitial,num
	FILEparam=open(sys.argv[2])
	linefileparam=FILEparam.readlines()
	dimparam=len(linefileparam)
	Jinitial=np.zeros(dimparam-1)	
	num=np.zeros(dimparam-1)
	i=0
	for line in linefileparam:
		if line.find('#') == -1:
			liner=line.split()
			Jinitial[i]=np.float32(liner[1])
			num[i]=np.float32(liner[0])
			print(num[i],Jinitial[i])
			i=i+1
	FILEparam.close()

def loadm():
	global m0,num2
	FILEparam2=open(sys.argv[2])
	linefileparam=FILEparam2.readlines()
	dimparam=len(linefileparam)
	m0=np.zeros(dimparam-1)	
	num2=np.zeros(dimparam-1)
	i=0
	for line in linefileparam:
		if line.find('#') == -1:
			liner=line.split()
			m0[i]=np.float32(liner[1])
			num2[i]=np.float32(liner[0])
			print(num2[i],m0[i])
			i=i+1
	FILEparam2.close()

def legendrefunc(x,A0,a2,a4):
	P0=legpoly(np.cos(x),[1,0,0,0,0])
	P2=legpoly(np.cos(x),[0,0,1,0,0])
	P4=legpoly(np.cos(x),[0,0,0,0,1])
	return A0*(P0+a2*P2+a4*P4)

#Get Coefficient of Legendre Fit from Experimental Intensity
def legendrefit():
	global popt,pcov,perr
	p0=np.array([float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5])])
	popt,pcov=cvt(legendrefunc,angle,intensity,p0)
	perr=np.sqrt(np.diag(pcov))
	print('A0:',popt[0],'a2:',popt[1],'a4:',popt[2])
	print('A0sderr:',perr[0],'a2sderr:',perr[1],'a4sderr:',perr[2])

#Experimental Intensity
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
def theo(x,delta,Ji,sigma):
	global Jf
	#Ji=int(sys.argv[5])
	Jf=float(sys.argv[6])
	#sigma=float(sys.argv[7])
	B0=ad.Bk(Ji,0,sigma)
	B2=ad.Bk(Ji,2,sigma)
	B4=ad.Bk(Ji,4,sigma)
	R00=ad.Rk(0,np.abs(Ji-Jf),np.abs(Ji-Jf),Ji,Jf)
	R01=ad.Rk(0,np.abs(Ji-Jf),np.abs(Ji-Jf)+1,Ji,Jf)
	R02=ad.Rk(0,np.abs(Ji-Jf)+1,np.abs(Ji-Jf)+1,Ji,Jf)
	R20=ad.Rk(2,np.abs(Ji-Jf),np.abs(Ji-Jf),Ji,Jf)
	R21=ad.Rk(2,np.abs(Ji-Jf),np.abs(Ji-Jf)+1,Ji,Jf)
	R22=ad.Rk(2,np.abs(Ji-Jf)+1,np.abs(Ji-Jf)+1,Ji,Jf)
	R40=ad.Rk(4,np.abs(Ji-Jf),np.abs(Ji-Jf),Ji,Jf)
	R41=ad.Rk(4,np.abs(Ji-Jf),np.abs(Ji-Jf)+1,Ji,Jf)
	R42=ad.Rk(4,np.abs(Ji-Jf)+1,np.abs(Ji-Jf)+1,Ji,Jf)
	Y0=B0*legpoly(np.cos(x),[1,0,0,0,0])*(R00+2.*R01*delta+R02*(delta**2.0))/(1.+(delta**2.))
	Y2=B2*legpoly(np.cos(x),[0,0,1,0,0])*(R20+2.*R21*delta+R22*(delta**2.0))/(1.+(delta**2.))
	Y4=B4*legpoly(np.cos(x),[0,0,0,0,1])*(R40+2.*R41*delta+R42*(delta**2.0))/(1.+(delta**2.))
	return Y0+Y2+Y4

#Normalized Theoretical Intensity w.r.t Experimental Intensity such that Yth(theta)=Yexp(theta):
def theonormmint(delta,Ji,sigma):
	numnorm=0.
	denumnorm=0.
	for l in range(0,dim-1,1):
		numnorm=numnorm+theo(angle[l],delta,Ji,sigma)*(intensity[l])/(error[l]**2.)
		denumnorm=denumnorm+((theo(angle[l],delta,Ji,sigma))**2.)/(error[l]**2.)
	#if round(denumnorm,3) > -0.001 and round(denumnorm,3) < 0.001:
		#ynorm=1.
	#else:
	ynorm=numnorm/denumnorm
	return ynorm

def chisq(delta,Ji,sigma):
	chi=0. #dim-2 is the number of data points -1
	for l in range(0,dim-1,1):
		chi=chi+(np.power(theonormmint(delta,Ji,sigma)*theo(angle[l],delta,Ji,sigma)-intensity[l],2.0))/((dim-2)*(error[l]**2.))
	return chi

def plotchi():
	#mixratio is angle for delta (in deg)
	mixratio=np.arange(-89,89,0.1)
	#delta in tan (degdel)
	delta=(np.tan(np.deg2rad(mixratio)))
	global mixratiomin, tandeltamin, chimin
	mixratiomin=np.zeros(3)
	tandeltamin=np.zeros(3)
	chimin=np.zeros(3)
	indchimin=np.zeros(3)
	fig2,ax2=plt.subplots()
	ax2.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax2.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax2.xaxis.set_minor_locator(tck.AutoMinorLocator())
	if option == 1:
		sigma=float(sys.argv[8])
		global msigma
		msigma = sigma
		for l in range(0,3,1):
			chimin[l]=np.min(chisq(delta,Jinitial[l],sigma))
			indchimin[l]=np.where(chisq(delta,Jinitial[l],sigma)==chimin[l])[0][0]
			mixratiomin[l]=mixratio[int(indchimin[l])]
			tandeltamin[l]=delta[int(indchimin[l])]
		ax2.plot(mixratio,chisq(delta,Jinitial[0],sigma),color='r',label='Ji:'+str(Jinitial[0])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'m:'+' '+str(round(sigma,3)))
		ax2.plot(mixratio,chisq(delta,Jinitial[1],sigma),color='b',label='Ji:'+str(Jinitial[1])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'m:'+' '+str(round(sigma,3)))
		ax2.plot(mixratio,chisq(delta,Jinitial[2],sigma),color='g',label='Ji:'+str(Jinitial[2])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'m:'+' '+str(round(sigma,3)))
		#print("Ji:", Jinitial[0],"delmin:",mixratio[np.where(chisq(delta,int(Jinitial[0]),sigma)==np.min(chisq(delta,int(Jinitial[0]),sigma)))[0][0]])
		#print("Ji:", Jinitial[1],"delmin:",mixratio[np.where(chisq(delta,int(Jinitial[1]),sigma)==np.min(chisq(delta,int(Jinitial[1]),sigma)))[0][0]])
		#print("Ji:", Jinitial[2],"delmin:",mixratio[np.where(chisq(delta,int(Jinitial[2]),sigma)==np.min(chisq(delta,int(Jinitial[2]),sigma)))[0][0]])
		print("Ji:", Jinitial[0],"delmin:", mixratiomin[0],"tandeltamin:", tandeltamin[0])
		print("Ji:", Jinitial[1],"delmin:", mixratiomin[1],"tandeltamin:", tandeltamin[1])
		print("Ji:", Jinitial[2],"delmin:", mixratiomin[2],"tandeltamin:", tandeltamin[2])
	else:
		Ji=float(sys.argv[8])
		global Jinit
		Jinit = Ji
		for l in range(0,3,1):
			chimin[l]=np.min(chisq(delta,Ji,m0[l]))
			indchimin[l]=np.where(chisq(delta,Ji,m0[l])==chimin[l])[0][0]
			mixratiomin[l]=mixratio[int(indchimin[l])]
			tandeltamin[l]=delta[int(indchimin[l])]
		ax2.plot(mixratio,chisq(delta,Ji,m0[0]),color='r',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'m:'+' '+str(round(m0[0],2)))
		ax2.plot(mixratio,chisq(delta,Ji,m0[1]),color='b',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'m:'+' '+str(round(m0[1],2)))
		ax2.plot(mixratio,chisq(delta,Ji,m0[2]),color='g',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'m:'+' '+str(round(m0[2],2)))
		#print("mi:", m0[0],"delmin:",mixratio[np.where(chisq(delta,Ji,m0[0])==np.min(chisq(delta,Ji,m0[0])))[0][0]])
		#print("mi:", m0[1],"delmin:",mixratio[np.where(chisq(delta,Ji,m0[1])==np.min(chisq(delta,Ji,m0[1])))[0][0]])
		#print("mi:", m0[2],"delmin:",mixratio[np.where(chisq(delta,Ji,m0[2])==np.min(chisq(delta,Ji,m0[2])))[0][0]])
		print("mi:", m0[0],"delmin:",mixratiomin[0],"tandeltamin:", tandeltamin[0])
		print("mi:", m0[1],"delmin:",mixratiomin[1],"tandeltamin:", tandeltamin[1])
		print("mi:", m0[2],"delmin:",mixratiomin[2],"tandeltamin:", tandeltamin[2])
	
	ax2.legend()
	#ax2.set_yscale("log")
	ax2.set_xlabel(r'$arctan(\delta)$',style='normal',fontweight='bold')
	ax2.set_ylabel('chisq',style='normal',fontweight='bold')
	ax2.set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	ax2.set_yscale('log')
	fig2.suptitle("32P Angular Distribution\nClarion2-Trinity\nO16+O18 at 30 MeV")

def writechiresult():
	with open(str(gamma)+'_ad_'+str(option)+'.txt','a') as fresults:
		fresults.write(str(dt.datetime.now())+'\n')
		if option == 1:
			fresults.write('m:'+' '+str(msigma)+'\n')
			fresults.write('#Ji\t#deltamin\t#tandeltamin\t#chi\n')
			for i in range(0,3,1):
				fresults.write(str(Jinitial[i])+'\t'+str(round(mixratiomin[i],4))+'\t'+str(round(tandeltamin[i],4))+'\t'+str(round(chimin[i],4))+'\n')
		else:
			fresults.write('Ji:'+' '+str(Jinit)+'\n')
			fresults.write('#mi\t#deltamin\t#tandeltamin\t#chi\n')
			for i in range(0,3,1):
				fresults.write(str(round(m0[i],2))+'\t'+str(round(mixratiomin[i],4))+'\t'+str(round(tandeltamin[i],4))+'\t'+str(round(chimin[i],4))+'\n')
		fresults.write('===========================================\n')

def main():
	loaddata()
	if option == 1:
		loadJi()
	else:
		loadm()
	legendrefit()
	plot()
	plotchi()
	writechiresult()
	plt.show()

if __name__ == "__main__":
	main()
