#!/usr/bin/env python3
#Legendre Fit for Angular Distribution as a Function of Angle
#Compare Theoretical Angular Distribution vs Mixing Ratio 
#C. Wibisono 
#06/29/'23
#07/11/'23 v2 with chisquare vs mixing ratio
#07/21/'24 Used this module for valid AD and Polarization Fit

#How to Use:
#./ad [addata] [list of Ji/list of sigma] A0_initial a2_initial a4_initial Jf [option] Ji/sigma [Qk factor]
#option == 1 for vary Ji fix m
#option == 2 for vary m fix Ji
#Note about Qk factor
#Qk factor == 1 for enabling qk factor into theoretical AD
#Qk factor == 0 for disabling qk factor

import numpy as np
from scipy.optimize import curve_fit as cvt
from numpy.polynomial.legendre import legval as legpoly
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import sys
import os
import adlib as ad
import adui
import polpredictlib2 as p
import datetime as dt
import time 
plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']
plt.rcParams['figure.dpi']= 200
#plt.rcParams['figure.autolayout']='True'
#plt.rcParams['figure.figsize']=(12,4)

option = int(sys.argv[7])
qkfactor = int(sys.argv[9])

def loaddata():
	global angle,intensity,error,gamma,dim,nuclei
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
			nuclei=liner2[3]
			print('Gamma Energy (keV):',gamma,'Nuclei:',nuclei)
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

def loaddetparam():
	global Radius,Distance,Thickness
	print("Enter Radius (cm):\n")
	Radius=float(input())
	print("Enter Distance (cm):\n")
	Distance=float(input())
	print("Enter Thickness (cm):\n")
	Thickness=float(input())

def loadpolparam():
	global thetapol,Jipol,Jfpol,mult,sigmapol,Pexp,Perr,polcolor,Jfsign,Plow,Pup
	print("Enter detector angle theta:\n")
	thetapol=float(input())
	print("Enter Ji:\n")
	Jipol=float(input())
	print("Enter Jf:\n")
	Jfpol=float(input())
	print("Enter sigma(magnetic width):\n")
	sigmapol=float(input())
	print("Enter polarization graph color: red (r), blue(b), green(g)\n")
	polcolor=input()
	print("Enter sign of Jf:\n")
	Jfsign=input()
	print("Enter L=Lo+1: 1 if E or 2 if M\n")
	mult=int(input())
	print("Enter Pmeasured\n")
	Pexp=float(input())
	print("Enter Perr\n")
	Perr=float(input())
	print("Enter lower limit of Polarization axis\n")
	Plow=float(input())
	print("Enter upper limit of Polarization axis\n")
	Pup=float(input())

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
	leg,=ax.plot(np.cos(xtrf),legendrefunc(xtrf,*popt),label='fit: A0=%5.3f, a2=%5.3f, a4=%5.3f' % tuple(popt))
	#ax.plot(np.cos(xtrf),legendrefunc(xtrf,11.09,-0.392,0.233),'m',label='AD')
	ax.legend()
	ax.set_xlabel(r'$cos(\theta)$',style='normal',fontweight='bold')
	ax.set_ylabel('Intensity(normalized)',style='normal',fontweight='bold')
	ax.set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	fig.suptitle(str(nuclei)+" "+"Angular Distribution\nClarion2-Trinity\nO16+O18 at 30 MeV")
	#fig.draw()
	line0=adui.adui(leg)
	line0.connect()
	plt.show()	


#Theoretical Intensity:
def theo(x,delta,Ji,sigma):
	global Jf
	#Ji=int(sys.argv[5])
	Jf=float(sys.argv[6])
	#sigma=float(sys.argv[7])
	B0=ad.Bk(Ji,0,sigma)
	B2=ad.Bk(Ji,2,sigma)
	B4=ad.Bk(Ji,4,sigma)
	if (abs(round(Ji,1)-round(Jf,1)) > 0) and (round(Jf,1) > 0):
		R00=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R01=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf)+1,Ji,Jf)
		R02=ad.Rk(0,abs(Ji-Jf)+1,abs(Ji-Jf)+1,Ji,Jf)
		R20=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R21=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf)+1,Ji,Jf)
		R22=ad.Rk(2,abs(Ji-Jf)+1,abs(Ji-Jf)+1,Ji,Jf)
		R40=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R41=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf)+1,Ji,Jf)
		R42=ad.Rk(4,abs(Ji-Jf)+1,abs(Ji-Jf)+1,Ji,Jf)
	
	elif (round(Ji,0) == 2 and round(Jf,0) == 0) or (round(Ji,0) == 1 and round(Jf,0) == 0) or (round(Ji,0) == 3 and round(Jf,0) == 0) or (round(Ji,1) == 0.5 and round(Jf,0)==0) or (round(Ji,1) == 1.5 and round(Jf,0) == 0) or (round(Ji,1) == 2.5 and round(Jf,0) == 0):
		R00=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R01=0
		R02=0
		R20=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R21=0
		R22=0
		R40=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R41=0
		R42=0
		
	else:
		R00=ad.Rk(0,1,1,Ji,Jf)
		R01=ad.Rk(0,1,2,Ji,Jf)
		R02=ad.Rk(0,2,2,Ji,Jf)
		R20=ad.Rk(2,1,1,Ji,Jf)
		R21=ad.Rk(2,1,2,Ji,Jf)
		R22=ad.Rk(2,2,2,Ji,Jf)
		R40=ad.Rk(4,1,1,Ji,Jf)
		R41=ad.Rk(4,1,2,Ji,Jf)
		R42=ad.Rk(4,2,2,Ji,Jf)
	
	if qkfactor == 1:
		Y0=B0*legpoly(np.cos(x),[1,0,0,0,0])*(R00+2.*R01*delta+R02*(delta**2.0))/(1.+(delta**2.))
		Y2=B2*ad.Qkcoeff(2,int(gamma),Radius,Distance,Thickness)*legpoly(np.cos(x),[0,0,1,0,0])*(R20+2.*R21*delta+R22*(delta**2.0))/(1.+(delta**2.))
		Y4=B4*ad.Qkcoeff(4,int(gamma),Radius,Distance,Thickness)*legpoly(np.cos(x),[0,0,0,0,1])*(R40+2.*R41*delta+R42*(delta**2.0))/(1.+(delta**2.))
	else:
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
	#if denumnorm > -0.0001 and denumnorm < 0.0001:
	#	ynorm=1.
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
	mixratio=np.arange(-89,89,1.) #used to overcome the entire mixing ratio (-inf,inf)-->equal to arctan(mixratio)
	#delta in tan (degdel)
	delta=(np.tan(np.deg2rad(mixratio))) #the true value of mixing ratio would be tan(arctan(mixratio))
	global mixratiomin, tandeltamin, chimin, chi1, chi2, chi3
	global fig2, ax2
	mixratiomin=np.zeros(3)
	tandeltamin=np.zeros(3)
	chimin=np.zeros(3)
	indchimin=np.zeros(3)
	#fig2,ax2=plt.subplots(2,1,figsize=(4.5,12),dpi=200,sharex=True)
	fig2,ax2=plt.subplots(2,1,figsize=(4.5,12),dpi=200,sharex=True)
	fig2.subplots_adjust(hspace=0)
	for j in range(0,2,1):
		if j == 0:
			ax2[j].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
			ax2[j].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
		else:
			ax2[j].tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='False',right='True',length=9,width=0.75)
			ax2[j].tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='False',right='True',length=6,width=0.75)
			ax2[j].tick_params(axis='y',labelsize='8')
		ax2[j].xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
		#ax2[j].yaxis.set_minor_locator(tck.AutoMinorLocator(n=5))

	if option == 1:
		sigma=float(sys.argv[8])
		global msigma
		msigma = sigma
		
		chi1,=ax2[0].plot(mixratio,chisq(delta,Jinitial[0],sigma),color='r',label='Ji:'+str(Jinitial[0])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+r'$\sigma$'+':'+' '+str(round(sigma,3)))
		chi2,=ax2[0].plot(mixratio,chisq(delta,Jinitial[1],sigma),color='b',label='Ji:'+str(Jinitial[1])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+r'$\sigma$'+':'+' '+str(round(sigma,3)))
		chi3,=ax2[0].plot(mixratio,chisq(delta,Jinitial[2],sigma),color='g',label='Ji:'+str(Jinitial[2])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+r'$\sigma$'+':'+' '+str(round(sigma,3)))
		if dim-1 == 4:
			ax2[0].axhline(y=5.423,linestyle='-.',color='k')
		if dim-1 == 5:
			ax2[0].axhline(y=4.616,linestyle='-.',color='k')
		if dim-1 == 6:
			ax2[0].axhline(y=4.103,linestyle='-.',color='k')
		if dim-1 == 7:
			ax2[0].axhline(y=3.743,linestyle='-.',color='k')
		if dim-1 == 8:
			ax2[0].axhline(y=3.475,linestyle='-.',color='k')
		if dim-1 == 9:
			ax2[0].axhline(y=3.266,linestyle='-.',color='k')
		#print("Ji:", Jinitial[0],"delmin:",mixratio[np.where(chisq(delta,int(Jinitial[0]),sigma)==np.min(chisq(delta,int(Jinitial[0]),sigma)))[0][0]])
		#print("Ji:", Jinitial[1],"delmin:",mixratio[np.where(chisq(delta,int(Jinitial[1]),sigma)==np.min(chisq(delta,int(Jinitial[1]),sigma)))[0][0]])
		#print("Ji:", Jinitial[2],"delmin:",mixratio[np.where(chisq(delta,int(Jinitial[2]),sigma)==np.min(chisq(delta,int(Jinitial[2]),sigma)))[0][0]])
		if round(Jf,0) > 0:
			for l in range(0,3,1):
				chimin[l]=np.min(chisq(delta,Jinitial[l],sigma))
				indchimin[l]=np.where(chisq(delta,Jinitial[l],sigma)==chimin[l])[0][0]
				mixratiomin[l]=mixratio[int(indchimin[l])]
				tandeltamin[l]=delta[int(indchimin[l])]
		
		print("Ji:", Jinitial[0],"delmin:", mixratiomin[0],"tandeltamin:", tandeltamin[0])
		print("Ji:", Jinitial[1],"delmin:", mixratiomin[1],"tandeltamin:", tandeltamin[1])
		print("Ji:", Jinitial[2],"delmin:", mixratiomin[2],"tandeltamin:", tandeltamin[2])
	else:
		Ji=float(sys.argv[8])
		global Jinit
		Jinit = Ji
		chi1,=ax2[0].plot(mixratio,chisq(delta,Ji,m0[0]),color='r',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+r'$\sigma$'+':'+' '+str(round(m0[0],2)))
		chi2,=ax2[0].plot(mixratio,chisq(delta,Ji,m0[1]),color='b',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+r'$\sigma$'+':'+' '+str(round(m0[1],2)))
		chi3,=ax2[0].plot(mixratio,chisq(delta,Ji,m0[2]),color='g',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+r'$\sigma$'+':'+' '+str(round(m0[2],2)))
		if dim-1 == 4:
			ax2[0].axhline(y=5.423,linestyle='-.',color='k')
		if dim-1 == 5:
			ax2[0].axhline(y=4.616,linestyle='-.',color='k')
		if dim-1 == 6:
			ax2[0].axhline(y=4.103,linestyle='-.',color='k')
		if dim-1 == 7:
			ax2[0].axhline(y=3.743,linestyle='-.',color='k')
		if dim-1 == 8:
			ax2[0].axhline(y=3.475,linestyle='-.',color='k')
		if dim-1 == 9:
			ax2[0].axhline(y=3.266,linestyle='-.',color='k')
		#print("mi:", m0[0],"delmin:",mixratio[np.where(chisq(delta,Ji,m0[0])==np.min(chisq(delta,Ji,m0[0])))[0][0]])
		#print("mi:", m0[1],"delmin:",mixratio[np.where(chisq(delta,Ji,m0[1])==np.min(chisq(delta,Ji,m0[1])))[0][0]])
		#print("mi:", m0[2],"delmin:",mixratio[np.where(chisq(delta,Ji,m0[2])==np.min(chisq(delta,Ji,m0[2])))[0][0]])
		if round(Jf,0) > 0:
			for l in range(0,3,1):
				chimin[l]=np.min(chisq(delta,Ji,m0[l]))
				indchimin[l]=np.where(chisq(delta,Ji,m0[l])==chimin[l])[0][0]
				mixratiomin[l]=mixratio[int(indchimin[l])]
				tandeltamin[l]=delta[int(indchimin[l])]
		print("sigma:", m0[0],"delmin:",mixratiomin[0],"tandeltamin:", tandeltamin[0])
		print("sigma:", m0[1],"delmin:",mixratiomin[1],"tandeltamin:", tandeltamin[1])
		print("sigma:", m0[2],"delmin:",mixratiomin[2],"tandeltamin:", tandeltamin[2])

	P=p.pol(thetapol,Jipol,Jfpol,mult,gamma,Pexp,Perr,sigmapol) #for Ji=5,6,7 choose the input of Ji as 6 so that the rest would be 5 and 7
	P2=p.pol(thetapol,Jipol,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
	
	if abs(round(Jipol,1)-round(Jfpol,1)) == 1 and Jfsign == '+':
		if mult == 1:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol+1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol+1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
			polplot1b,=ax2[1].plot(mixratio,Pb,color='r',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='g',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='g',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
		if mult == 2:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol+1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol+1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot1b,=ax2[1].plot(mixratio,Pb,color='r',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='g',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='g',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
	
	if abs(round(Jipol,1)-round(Jfpol,1)) == 1 and Jfsign == '-':
		if mult == 1:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol+1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol+1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
		
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot1b,=ax2[1].plot(mixratio,Pb,color='r',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='g',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='g',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
		if mult == 2:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol+1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol+1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot1b,=ax2[1].plot(mixratio,Pb,color='r',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol-1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='g',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='g',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol+1)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
	
	if abs(round(Jipol,1)-round(Jfpol,1)) == 2 and Jfsign == '+':
		if mult == 1:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol-2,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol-2,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M2-E3
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E2-M3
			polplot1b,=ax2[1].plot(mixratio,Pb,color='b',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='b',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='r',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2

		if mult == 2:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol-2,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol-2,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E2-M3
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M2-E3
			polplot1b,=ax2[1].plot(mixratio,Pb,color='b',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='b',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='r',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E1-M2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M1-E2

	if abs(round(Jipol,1)-round(Jfpol,1)) == 2 and Jfsign == '-':
		
		if mult == 1:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol-2,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol-2,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M2-E3
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E2-M3
			polplot1b,=ax2[1].plot(mixratio,Pb,color='b',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='b',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='r',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
		
		if mult == 2:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol-2,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol-2,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E2-M3
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M2-E3
			polplot1b,=ax2[1].plot(mixratio,Pb,color='b',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot2b,=ax2[1].plot(mixratio,P2b,color='b',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
			polplot1c,=ax2[1].plot(mixratio,Pc,color='r',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E1-M2
			polplot2c,=ax2[1].plot(mixratio,P2c,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M1-E2
	
	if abs(round(Jipol,1)-round(Jfpol,1)) == 3 and Jfsign == '+':
		if mult == 1:
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M3-E4
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E3-M4
		if mult == 2:
			Pb=p.pol(thetapol,Jipol-1,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2b=p.pol(thetapol,Jipol-1,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			Pc=p.pol(thetapol,Jipol-2,Jfpol,mult,gamma,Pexp,Perr,sigmapol)
			P2c=p.pol(thetapol,Jipol-2,Jfpol,mult-1,gamma,Pexp,Perr,sigmapol)
			
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E2-M3
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M2-E3
			polplot1b,=ax2[1].plot(mixratio,Pb,color='b',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E2-M3
			polplot2b,=ax2[1].plot(mixratio,P2b,color='b',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M2-E3
			polplot1c,=ax2[1].plot(mixratio,Pc,color='r',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #E2-M3
			polplot2c,=ax2[1].plot(mixratio,P2c,color='r',linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{+}$') #M2-E3

	


	if abs(round(Jipol,1)-round(Jfpol,1)) == 3 and Jfsign == '-':
		if mult == 1:
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M3-E4
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E3-M4
		if mult == 2:
			polplot1,=ax2[1].plot(mixratio,P,color=polcolor,label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{+}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #E3-M4
			polplot2,=ax2[1].plot(mixratio,P2,color=polcolor,linestyle='-.',label='theory:'+' '+'Ji:'+str(Jipol)+r'$^{-}$'+'----->'+'Jf:'+str(Jfpol)+r'$^{-}$') #M3-E4

	
	ax2[1].axhline(y=Pexp,color='k',label='experiment')
	ax2[1].axhspan(Pexp-Perr,Pexp+Perr,alpha=0.25,color='g')

	#ax2[0].legend(loc='upper right')
	#ax2[1].legend(loc='upper right')
	#ax2.set_yscale("log")
	ax2[1].set_xlabel(r'$tan^{-1} (\delta)$'+' '+'(degrees)',style='normal',fontweight='bold')
	ax2[0].set_ylabel(r'$\chi^{2}$',style='normal',fontweight='bold')
	ax2[1].set_ylabel('Polarization',style='normal',fontsize='small',fontweight='bold')
	#ax2[0].set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	ax2[0].set_yscale('log')
	ax2[1].set_ylim(Plow,Pup)
	ax2[1].set_xlim(-89,89)
	ax2[0].set_xlim(-89,89)
	#fig2.suptitle(str(nuclei)+" "+"Angular Distribution and Polarization Measurement\n"+"Gamma Energy:"+" "+str(gamma)+" "+"keV"+"\n"+"Clarion2-Trinity: O16+O18 at 30 MeV")
	#plt.savefig('832file.eps',dpi='figure',format='eps')
	#fig2.canvas.draw()
	line1=adui.adui(chi1)
	line2=adui.adui(chi2)
	line3=adui.adui(chi3)
	line10=adui.adui(polplot1)
	line1.connect()
	line2.connect()
	line3.connect()
	line10.connect()
	plt.show()


#Delta to plot theoretical intensity:
def loaddelta():
	global deltaoption,tandeltamod
	tandeltamod = np.zeros(3)
	print("Use global deltamin: Yes(1) No(0)\n")
	deltaoption = int(input())
	if deltaoption == 0:
		print("Manually Enter delta(?): Yes(1) No(0)\n")
		enterdelta=int(input())
		if enterdelta == 1:
			print("Input arc tan delta min\n")
			print("Enter delta for Ji1/sigmai1\n")
			tandeltamod[0] = np.tan(np.deg2rad(float(input())))
			print("Enter delta for Ji2/sigmai2\n")
			tandeltamod[1] = np.tan(np.deg2rad(float(input())))
			print("Enter delta for Ji3/sigmai3\n")
			tandeltamod[2] = np.tan(np.deg2rad(float(input())))
		else:
			print("Rerun the program to generate chi-square and record the local minima for each plot\n")
			exit()



#Plot theoretical Angular Distribution Fit along with Data	
def plotad():
	fig3,ax3=plt.subplots()
	ax3.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
	ax3.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
	ax3.xaxis.set_minor_locator(tck.AutoMinorLocator(n=5))
	angledeg=np.arange(0,181,0.01)
	anglerad=np.deg2rad(angledeg)
	if option == 1:
		sigma=float(sys.argv[8])

		if deltaoption == 1: #plot theoAD with global minimum of delta
			ax3.plot(angledeg,theonormmint(tandeltamin[0],Jinitial[0],sigma)*theo(anglerad,tandeltamin[0],Jinitial[0],sigma),color='r',label='Ji:'+str(Jinitial[0])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'delta:'+' '+str(round(tandeltamin[0],3)))
			ax3.plot(angledeg,theonormmint(tandeltamin[1],Jinitial[1],sigma)*theo(anglerad,tandeltamin[1],Jinitial[1],sigma),color='b',label='Ji:'+str(Jinitial[1])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'delta:'+' '+str(round(tandeltamin[1],3)))
			ax3.plot(angledeg,theonormmint(tandeltamin[2],Jinitial[2],sigma)*theo(anglerad,tandeltamin[2],Jinitial[2],sigma),color='g',label='Ji:'+str(Jinitial[2])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'delta:'+' '+str(round(tandeltamin[2],3)))

		if deltaoption == 0: #plot theoAD based on user's delta
			ax3.plot(angledeg,theonormmint(tandeltamod[0],Jinitial[0],sigma)*theo(anglerad,tandeltamod[0],Jinitial[0],sigma),color='r',label='Ji:'+str(Jinitial[0])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'delta:'+' '+str(round(tandeltamod[0],3)))
			ax3.plot(angledeg,theonormmint(tandeltamod[1],Jinitial[1],sigma)*theo(anglerad,tandeltamod[1],Jinitial[1],sigma),color='b',label='Ji:'+str(Jinitial[1])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'delta:'+' '+str(round(tandeltamod[1],3)))
			ax3.plot(angledeg,theonormmint(tandeltamod[2],Jinitial[2],sigma)*theo(anglerad,tandeltamod[2],Jinitial[2],sigma),color='g',label='Ji:'+str(Jinitial[2])+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'delta:'+' '+str(round(tandeltamod[2],3)))
	
		ax3.errorbar(np.rad2deg(angle),intensity,yerr=error,fmt='o',linewidth=2,capsize=6,label='data')

	else:
		Ji=float(sys.argv[8])
		
		if deltaoption == 1: #plot theoAD based on global minimum of delta
			ax3.plot(angledeg,theonormmint(tandeltamin[0],Ji,m0[0])*theo(anglerad,tandeltamin[0],Ji,m0[0]),color='r',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'sigma:'+' '+str(round(m0[0],2))+' '+'delta:'+str(round(tandeltamin[0],3)))
			ax3.plot(angledeg,theonormmint(tandeltamin[1],Ji,m0[1])*theo(anglerad,tandeltamin[1],Ji,m0[1]),color='b',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'sigma:'+' '+str(round(m0[1],2))+' '+'delta:'+str(round(tandeltamin[1],3)))
			ax3.plot(angledeg,theonormmint(tandeltamin[2],Ji,m0[2])*theo(anglerad,tandeltamin[2],Ji,m0[2]),color='g',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'sigma:'+' '+str(round(m0[2],2))+' '+'delta:'+str(round(tandeltamin[2],3)))
		
		if deltaoption == 0: #plot theoAD based on user's delta
			
			ax3.plot(angledeg,theonormmint(tandeltamod[0],Ji,m0[0])*theo(anglerad,tandeltamod[0],Ji,m0[0]),color='r',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'sigma:'+' '+str(round(m0[0],2))+' '+'delta:'+str(round(tandeltamod[0],3)))
			ax3.plot(angledeg,theonormmint(tandeltamod[1],Ji,m0[1])*theo(anglerad,tandeltamod[1],Ji,m0[1]),color='b',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'sigma:'+' '+str(round(m0[1],2))+' '+'delta:'+str(round(tandeltamod[1],3)))
			ax3.plot(angledeg,theonormmint(tandeltamod[2],Ji,m0[2])*theo(anglerad,tandeltamod[2],Ji,m0[2]),color='g',label='Ji:'+str(Ji)+' '+'--->'+' '+'Jf:'+str(Jf)+' '+'sigma:'+' '+str(round(m0[2],2))+' '+'delta:'+str(round(tandeltamod[2],3)))
		
		ax3.errorbar(np.rad2deg(angle),intensity,yerr=error,fmt='o',linewidth=2,capsize=6,label='data')

	ax3.legend()
	#ax2.set_yscale("log")
	ax3.set_xlabel(r'$\theta_{det} (deg)$',style='normal',fontweight='bold')
	ax3.set_ylabel('Intensity',style='normal',fontweight='bold')
	ax3.set_title('Gamma Energy:'+' '+str(gamma)+' '+'keV')
	fig3.suptitle(str(nuclei)+" "+"Angular Distribution\nClarion2-Trinity\nO16+O18 at 30 MeV")
	#plt.savefig('1555file.eps',dpi='figure',format='eps')
	#fig3.canvas.draw()

def theoa2a4(Ji,delta,sigma):
	B0=ad.Bk(Ji,0,sigma)
	B2=ad.Bk(Ji,2,sigma)
	B4=ad.Bk(Ji,4,sigma)
	if abs(round(Ji,1)-round(Jf,1)) > 0:
		R00=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R01=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf)+1,Ji,Jf)
		R02=ad.Rk(0,abs(Ji-Jf)+1,abs(Ji-Jf)+1,Ji,Jf)
		R20=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R21=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf)+1,Ji,Jf)
		R22=ad.Rk(2,abs(Ji-Jf)+1,abs(Ji-Jf)+1,Ji,Jf)
		R40=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R41=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf)+1,Ji,Jf)
		R42=ad.Rk(4,abs(Ji-Jf)+1,abs(Ji-Jf)+1,Ji,Jf)
	
	elif round(Ji,0) == 2 and round(Jf,0) == 0:
		R00=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R01=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R02=ad.Rk(0,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R20=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R21=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R22=ad.Rk(2,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R40=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R41=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		R42=ad.Rk(4,abs(Ji-Jf),abs(Ji-Jf),Ji,Jf)
		
	else:
		R00=ad.Rk(0,1,1,Ji,Jf)
		R01=ad.Rk(0,1,2,Ji,Jf)
		R02=ad.Rk(0,2,2,Ji,Jf)
		R20=ad.Rk(2,1,1,Ji,Jf)
		R21=ad.Rk(2,1,2,Ji,Jf)
		R22=ad.Rk(2,2,2,Ji,Jf)
		R40=ad.Rk(4,1,1,Ji,Jf)
		R41=ad.Rk(4,1,2,Ji,Jf)
		R42=ad.Rk(4,2,2,Ji,Jf)
	
	if qkfactor == 1:
		a0th=B0*(R00+2.*R01*delta+R02*(delta**2.0))/(1.+(delta**2.))
		a2th=B2*ad.Qkcoeff(2,int(gamma),Radius,Distance,Thickness)*(R20+2.*R21*delta+R22*(delta**2.0))/(1.+(delta**2.))
		a4th=B4*ad.Qkcoeff(4,int(gamma),Radius,Distance,Thickness)*(R40+2.*R41*delta+R42*(delta**2.0))/(1.+(delta**2.))

	else:
		a0th=B0*(R00+2.*R01*delta+R02*(delta**2.0))/(1.+(delta**2.))
		a2th=B2*(R20+2.*R21*delta+R22*(delta**2.0))/(1.+(delta**2.))
		a4th=B4*(R40+2.*R41*delta+R42*(delta**2.0))/(1.+(delta**2.))
	
	return a0th,a2th,a4th

def evaltheoa2a4():
	global a0theo,a2theo,a4theo
	a0theo=np.zeros(3)
	a2theo=np.zeros(3)
	a4theo=np.zeros(3)
	if option == 1:
		if deltaoption == 1: #plot theoAD based on global minimum of delta:
			for i in range(0,3,1):
				a0theo[i],a2theo[i],a4theo[i]=theoa2a4(Jinitial[i],tandeltamin[i],msigma)
		if deltaoption == 0:
			for i in range(0,3,1):
				a0theo[i],a2theo[i],a4theo[i]=theoa2a4(Jinitial[i],tandeltamod[i],msigma)
	
		print("Theoretical a0,a2,a4:\n")
		for i in range(0,3,1):
			print("Ji:", Jinitial[i],"a0th:", a0theo[i],"a2th:", a2theo[i],"a4th:", a4theo[i]) 
	else:
		if deltaoption == 1:
			for i in range(0,3,1):
				a0theo[i],a2theo[i],a4theo[i]=theoa2a4(Jinit,tandeltamin[i],m0[i])		
		if deltaoption == 0:
			for i in range(0,3,1):
				a0theo[i],a2theo[i],a4theo[i]=theoa2a4(Jinit,tandeltamod[i],m0[i])
		print("Theoretical a0,a2,a4:\n")
		for i in range(0,3,1):
			print("sigma:", m0[i],"a0th:", a0theo[i],"a2th:", a2theo[i],"a4th:", a4theo[i]) 


def writechiresult():
	with open(str(gamma)+'_ad_'+str(option)+'.txt','a') as fresults:
		fresults.write(str(dt.datetime.now())+'\n')

		if qkfactor == 1:
			fresults.write('Attenuation Factor:'+' '+ 'yes\n')
			fresults.write('Radius:'+' '+str(Radius)+' '+'cm\n')
			fresults.write('Distance:'+' '+str(Distance)+' '+'cm\n')
			fresults.write('Thickness:'+' '+str(Thickness)+' '+'cm\n')
		else:
			fresults.write('Attenuation Factor:'+' '+'no\n')

		fresults.write('Data:\n')
		fresults.write('Gamma Energy:'+' '+str(gamma)+'\n')
		fresults.write('Angle(Deg)'+'\t'+'Intensity'+'\t'+'Error'+'\n')
		for i in range(0,dim-1,1):
			fresults.write(str(round(np.rad2deg(angle[i]),2))+'\t'+str(round(intensity[i],3))+'\t'+str(round(error[i],3))+'\n')
		fresults.write('Legendre Fit:\n')
		fresults.write('A0\ta2\ta4\n')
		fresults.write(str(round(popt[0],4))+'\t'+str(round(popt[1],4))+'\t'+str(round(popt[2],4))+'\n')
		fresults.write('stdA0\tstda2\tstda4\n')
		fresults.write(str(round(perr[0],4))+'\t'+str(round(perr[1],4))+'\t'+str(round(perr[2],4))+'\n')

		if option == 1:
			fresults.write('sigma:'+' '+str(msigma)+'\n')
			fresults.write('Theoretical a0,a2,a4:\n')
			fresults.write('Ji\ta0th\ta2th\ta4th\n')
			for i in range(0,3,1):
				fresults.write(str(Jinitial[i])+'\t'+str(round(a0theo[i],4))+'\t'+str(round(a2theo[i],4))+'\t'+str(round(a4theo[i],4))+'\n')
			if deltaoption == 1:
				fresults.write('#Ji\t#arctan(deltamin)\t#tanarctan(deltamin)\t#chi\n')
				for i in range(0,3,1):
					fresults.write(str(Jinitial[i])+'\t'+str(round(mixratiomin[i],4))+'\t'+str(round(tandeltamin[i],4))+'\t'+str(round(chimin[i],4))+'\n')
			if deltaoption == 0:
				fresults.write('#Ji\t#arctan(deltamin)\t#tanarctan(deltamin)\n')
				for i in range(0,3,1):
					fresults.write(str(Jinitial[i])+'\t'+str(round(np.rad2deg(np.arctan(tandeltamod[i])),4))+'\t'+str(round(tandeltamod[i],4))+'\n')
				
		else:
			fresults.write('Ji:'+' '+str(Jinit)+'\n')
			fresults.write('Theoretical a0,a2,a4:\n')
			fresults.write('mi\ta0th\ta2th\ta4th\n')
			for i in range(0,3,1):
				fresults.write(str(round(m0[i],2))+'\t'+str(round(a0theo[i],4))+'\t'+str(round(a2theo[i],4))+'\t'+str(round(a4theo[i],4))+'\n')
			if deltaoption == 1:
				fresults.write('#sigmai\t#arctan(deltamin)\t#tanarctan(deltamin)\t#chi\n')
				for i in range(0,3,1):
					fresults.write(str(round(m0[i],2))+'\t'+str(round(mixratiomin[i],4))+'\t'+str(round(tandeltamin[i],4))+'\t'+str(round(chimin[i],4))+'\n')
			if deltaoption == 0:
				fresults.write('#sigmai\t#arctan(deltamin)\t#tanarctan(deltamin)\n')
				for i in range(0,3,1):
					fresults.write(str(round(m0[i],2))+'\t'+str(round(np.rad2deg(np.arctan(tandeltamod[i])),4))+'\t'+str(round(tandeltamod[i],4))+'\n')
		fresults.write('===========================================\n')

def main():
	start=time.time()
	loaddata()
	if option == 1:
		loadJi()
	else:
		loadm()
	if qkfactor == 1:
		loaddetparam()

	legendrefit()
	plot()
	loadpolparam()
	print("Generating Chi-Square Plot....\n")
	print("Generating Polarization Plot....\n")
	plotchi()
	#loaddelta()
	#print("Generating Theoritical AD Plot...\n")
	#plotad()
	#print("Calculating Theoretical a2 and a4...\n")
	#evaltheoa2a4()
	writechiresult()
	end=time.time()
	print("The execution time:", (end-start), "s")
	plt.show()

if __name__ == "__main__":
	main()
