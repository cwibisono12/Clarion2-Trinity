#!/usr/bin/env python3

#./vis.py file 5000 5000 0 4999 0 4999 
#C. Wibisono

#Program for Displaying 2D Histogram + Drawing Banana Cuts + Making Basic Projections:

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import sys
import os
from matplotlib.colors import LogNorm
from scipy.optimize import curve_fit as cvt


print("Welcome to Clarion2-Trinity Portable Visualization Software from C.W\n")
print("Here are the features of the Software:")
print("Hit g to make Projection")
print("Hit e to expand axis")
print("Hit f to generate coordinates useful for making banana gates")
print("click mouse to draw banana gates")
print("Hit n to perform Gauss fit")
print("Hit m to perform another type of Gauss fit")
print("Hit d to perform Double Gaussian")
print("Hit q to close figure")

#Reading Matrix Files:
y=np.fromfile(sys.argv[1],dtype=np.int32,sep="",count=-1)

#Transpose 1D to 2D:
ytrp=np.reshape(y,(int(sys.argv[2]),int(sys.argv[3])))

#2D-Matrix:
fig, ax=plt.subplots()

pos=ax.imshow(ytrp,cmap="gist_ncar",origin="lower",extent=[int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7])],aspect="auto",norm=LogNorm())

plt.colorbar(pos,ax=ax)


#1D-Matrix:
figp, axp=plt.subplots()

axp.tick_params(direction='in',axis='both',which='major',bottom='True',left='True',top='True',right='True',length=9,width=0.75)
axp.tick_params(direction='in',axis='both',which='minor',bottom='True',left='True',top='True',right='True',length=6,width=0.75)
plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

#axp=[]
#Prepared the array for y and x:
yproj=np.zeros(int(sys.argv[2]),dtype=np.int32) #Project onto Y axis
xproj=np.zeros(int(sys.argv[3]),dtype=np.int32) #Project onto X axis

yprojb=np.zeros(int(sys.argv[2]),dtype=np.int32) #yarr for background
xprojb=np.zeros(int(sys.argv[3]),dtype=np.int32) #xarr for background

yprojbs=np.zeros(int(sys.argv[2]),dtype=np.int32) #yarr for background subtracted
xprojbs=np.zeros(int(sys.argv[3]),dtype=np.int32) #xarr for background subtracted

xx=np.zeros(int(sys.argv[3]),dtype=np.int32)
xy=np.zeros(int(sys.argv[2]),dtype=np.int32)

#Gaussian Function:
def gauss(x,H,A,mu,sigma):
	return H+A*np.exp(-((x-mu)**2.)/(2*(sigma**2.)))

backgheight=0
def gauss2(x,A,mu,sigma,backgheight):
	return A*np.exp(-((x-mu)**2.)/(2*(sigma**2.)))+backgheight
	
def gausspure(x,A,mu,sigma):
	return A*np.exp(-((x-mu)**2.)/(2*(sigma**2.)))

#Double Gaussian Function:
#def gaussdoub(x,H1,A1,mu1,sigma1,H2,A2,mu2,sigma2):
#	return gauss(x,H1,A1,mu1,sigma1)+gauss(x,H2,A2,mu2,sigma2)

backgheight2=0
#Type2:
def gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight2):
	return gausspure(x,A1,mu1,sigma1)+gausspure(x,A2,mu2,sigma2)+backgheight2

coordsx=[]
coordsy=[]

index=0
proj=0


def onclick(event):
	global index
	index=index+1
	xcoords, ycoords=event.xdata, event.ydata
	xpix, ypix=event.x,event.y
	global coords
	coordsx.append(xcoords)
	coordsy.append(ycoords)
	ax.plot(coordsx,coordsy,'ro-')
	fig.canvas.draw()
	print(index,xcoords,ycoords,xpix,ypix)

def onpress(event):
	##Zoom in and Zoom out Histogram:
	if event.key == 'e':
		print("Enter Axes to Expand 1 or 2\n")
		waxis=input()
		global axp,xx,xy,proj,background,xlow,xup
		global xproj,yproj,xprojbs,yprojbs
		global low,high,gate
		if int(waxis) == 1:
			print("Enter the lower xlim\n")
			xlow1=input()
			print("Enter the upper xlim\n")
			xup1=input()
			print("Enter the lower ylim\n")
			ylow1=input()
			print("Enter the upper ylim\n")
			yup1=input()
			ax.set_xlim(int(xlow1),int(xup1))
			ax.set_ylim(int(ylow1),int(yup1))
		if int(waxis) == 2: 
			print("Enter the lower limit\n")
			xlow=input()
			print("Enter the upper limit\n")
			xup=input()
			axp.clear()
			axp.set_xlim(int(xlow),int(xup))
			if int(proj)==1 and int(background)==0:	
				axp.plot(xx[int(xlow):int(xup)],xproj[int(xlow):int(xup)],linewidth=0.75,ls='steps',label='ProjX')
				axp.legend()
				axp.set_title('Full Projection X')
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
			if int(proj)==1 and int(background)==1:	
				axp.plot(xx[int(xlow):int(xup)],xprojbs[int(xlow):int(xup)],linewidth=0.75,ls='steps',label='ProjX')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Gate on '+str(gate)+' keV')
			if int(proj)==2 and int(background)==0:
				axp.plot(xy[int(xlow):int(xup)],yproj[int(xlow):int(xup)],linewidth=0.75,ls='steps',label='ProjY')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection Y')
			if int(proj)==2 and int(background)==1:
				axp.plot(xy[int(xlow):int(xup)],yprojbs[int(xlow):int(xup)],linewidth=0.75,ls='steps',label='ProjY')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Gate on '+str(gate)+' keV')
			plt.draw()
	#Zoom out Histogram:
	if event.key == 'o':
		ax.set_xlim(0,int(sys.argv[3]))
		ax.set_ylim(0,int(sys.argv[2]))
		if int(proj)==1:
			axp.set_xlim(0,int(sys.argv[3]))
		if int(proj)==2:
			axp.set_xlim(0,int(sys.argv[2]))
	
	#Saving Coordinates File for Banana Gates:
	if event.key == 'f':
		print("Print Ban Coords to the Screen:\n")
		print("Input Banana ID:\n")
		banid=input()
		print(banid,index,"#")
		zx=np.zeros(index)
		zy=np.zeros(index)
		for i in range(index):
			zx[i]=coordsx[i]
			zy[i]=coordsy[i]
			print(banid,index,zx[i],zy[i])
	
	if event.key == 'g':
		print("ProjectionX or ProjectionY 1 or 2:\n")
		proj=input()
		print("Enter the region where you want to make gates:\n")
		print("Enter the lower region:\n")
		low=input()
		print("Enter the higher region:\n")
		high=input()
		widthgate=int(high)-int(low)
		gate=int(int(low)+widthgate/2.)
		#global axp
		#Projection on X axis: Make a gate on Y axis then project on X axis:
		axp.clear()
		if int(proj) == 1:
			#global x
			xx=np.arange(0,int(sys.argv[3]),1)
			xproj=np.zeros(int(sys.argv[3]),dtype=np.int32)
			xprojb=np.zeros(int(sys.argv[3]),dtype=np.int32)
			xprojbs=np.zeros(int(sys.argv[3]),dtype=np.int32)
			print("Perform Background Subtraction Yes or No hit 0 or 1:\n")
			background=input()
			
			#Create X Projection from Y Gate:
			for i in range(0,int(sys.argv[3]),1):
				for j in range(int(low),int(high)+1,1):
					xproj[i]=xproj[i]+ytrp[j,i]
			if int(background) == 0:	
				axp.plot(xx,xproj,linewidth=0.75,ls='steps',label='ProjX')
				axp.set_title('Full Projection X')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
			if int(background) == 1:
				print("Enter the left region for background")
				backleft=input()
				print("Enter the right region for background")
				backright=input()
				widthbackg=int(backright)-int(backleft)
				for i in range(0,int(sys.argv[2]),1):
					for j in range(int(backleft),int(backright)+1,1):
						xprojb[i]=xprojb[i]+ytrp[j,i]
					xprojbs[i]=(widthgate/(widthbackg+widthgate))*xproj[i]-(widthbackg/(widthbackg+widthgate))*xprojb[i]
				axp.clear()
				axp.plot(xx,xprojbs,linewidth=0.75,ls='steps',label='ProjX')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())


		#Projection on Y axis: Make a gate on X axis then project on Y axis:
		if int(proj) == 2:
			#global x
			xy=np.arange(0,int(sys.argv[2]),1)
			yproj=np.zeros(int(sys.argv[2]),dtype=np.int32)
			yprojb=np.zeros(int(sys.argv[2]),dtype=np.int32)
			yprojbs=np.zeros(int(sys.argv[2]),dtype=np.int32)
			print("Perform Background Subtraction Yes or No hit 0 or 1:\n")
			background=input()
			#Create Y Projection from X Gate:
			for j in range(0,int(sys.argv[2]),1):
				for i in range(int(low),int(high)+1,1):
					yproj[j]=yproj[j]+ytrp[j,i]
			if int(background) == 0:
				axp.plot(xy,yproj,linewidth=0.75,ls='steps',label='ProjY')
				axp.set_title('Full Projection Y')
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.legend()
			if int(background) == 1:
				print("Enter the left region for background")
				backleft=input()
				print("Enter the right region for background")
				backright=input()
				widthbackg=int(backright)-int(backleft)
				for j in range(0,int(sys.argv[2]),1):
					for i in range(int(backleft),int(backright)+1,1):
						yprojb[j]=yprojb[j]+ytrp[j,i]
					yprojbs[j]=(widthgate*yproj[j]-widthbackg*yprojb[j])/(widthgate+widthbackg)
				axp.clear()
				axp.plot(xy,yprojbs,linewidth=0.75,ls='steps',label='ProjY')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
	
	#Perform Gauss fit for a particular peak:
	if event.key == 'n':
		print("Perform Gauss Fit:\n")
		print("Set the background height:\n")
		backgfit=input()
		print("Set the approximate peak height:\n")
		height=input()
		print("Set the approximate center position:\n")
		center=input()
		print("Set the approximate width:\n")
		std=input()
		print("Select the region of interest:\n")
		print("select lower xlim:\n")
		xlowg=input()
		print("select the upper xlim:\n")
		xupg=input()
		p0=np.array([int(backgfit),int(height),int(center),int(std)])
		axp.clear()
		if int(proj) == 1:
			if int(background) == 0:
				axp.plot(xx[int(xlow):int(xup)+1],xproj[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjX')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection X')
				popt,pcov=cvt(gauss,xx[int(xlowg):int(xupg)+1],xproj[int(xlowg):int(xupg)+1],p0)

			if int(background) == 1:	
				axp.plot(xx[int(xlow):int(xup)+1],xprojbs[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjX')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				popt,pcov=cvt(gauss,xx[int(xlowg):int(xupg)+1],xprojbs[int(xlowg):int(xupg)+1],p0)
	
			axp.plot(xx[int(xlowg):int(xupg)+1],gauss(xx[int(xlowg):int(xupg)+1],*popt),'r',linewidth=0.5)
			area=(np.sqrt(2*(np.pi)))*popt[1]*abs(popt[3])
			print("mean:",popt[2],"sigma:",abs(popt[3]),"area:",area)

		if int(proj) == 2:
			if int(background) == 0:
				axp.plot(xy[int(xlow):int(xup)+1],yproj[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjY')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection Y')
				popt,pcov=cvt(gauss,xy[int(xlowg):int(xupg)+1],yproj[int(xlowg):int(xupg)+1],p0)

			if int(background) == 1:	
				axp.plot(xy[int(xlow):int(xup)+1],yprojbs[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjY')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				popt,pcov=cvt(gauss,xy[int(xlowg):int(xupg)+1],yprojbs[int(xlowg):int(xupg)+1],p0)	
			
			axp.plot(xy[int(xlowg):int(xupg)+1],gauss(xy[int(xlowg):int(xupg)+1],*popt),'r',linewidth=0.5)
			area=(np.sqrt(2*(np.pi)))*popt[1]*abs(popt[3])
			print("mean:",popt[2],"sigma:",abs(popt[3]),"area:",area)
	
#Perform another type of Gauss fit:
	if event.key == 'm':
		print("Perform Gauss Fit:\n")
		print("Set the background left:\n")
		backgfitleft=input()
		print("Set the background right:\n")
		backgfitright=input()
		print("Select the region of interest:\n")
		print("select lower xlim:\n")
		xlowg=input()
		print("select the upper xlim:\n")
		xupg=input()
		backgfitcenter=int(backgfitleft)+int((int(backgfitright)-int(backgfitleft))/2)
		meanfitcenter=int(xlowg)+int((int(xupg)-int(xlowg))/2)
		axp.clear()
		if int(proj) == 1:
			if int(background) == 0:
				axp.plot(xx[int(xlow):int(xup)+1],xproj[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjX')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection X')
				p0=np.array([xproj[int(meanfitcenter)],xx[int(meanfitcenter)],xx[2*(int(meanfitcenter)-int(xlowg))]])
				backgheight=np.average(xproj[int(backgfitleft):int(backgfitright)])
				print(backgfitcenter,backgheight)
				popt,pcov=cvt(lambda x, A, mu, sigma: gauss2(x,A,mu,sigma,backgheight),xx[int(xlowg):int(xupg)+1],xproj[int(xlowg):int(xupg)+1],p0)

			if int(background) == 1:	
				axp.plot(xx[int(xlow):int(xup)+1],xprojbs[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjX')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				p0=np.array([xprojbs[int(meanfitcenter)],xx[int(meanfitcenter)],xx[2*(int(meanfitcenter)-int(xlowg))]])
				backgheight=np.average(xprojbs[int(backgfitleft):int(backgfitright)])
				print(backgfitcenter,backgheight)
				popt,pcov=cvt(lambda x, A, mu, sigma: gauss2(x,A,mu,sigma,backgheight),xx[int(xlowg):int(xupg)+1],xprojbs[int(xlowg):int(xupg)+1],p0)
	
			axp.plot(xx[int(xlowg):int(xupg)+1],gauss2(xx[int(xlowg):int(xupg)+1],*popt,backgheight),'r',linewidth=0.5)
			area=(np.sqrt(2*(np.pi)))*popt[0]*abs(popt[2])
			print("mean:",popt[1],"sigma:",abs(popt[2]),"area:",area)

		if int(proj) == 2:
			if int(background) == 0:
				axp.plot(xy[int(xlow):int(xup)+1],yproj[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjY')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection Y')
				p0=np.array([yproj[int(meanfitcenter)],xy[int(meanfitcenter)],xy[2*(int(meanfitcenter)-int(xlowg))]])
				backgheight=np.average(yproj[int(backgfitleft):int(backgfitright)])
				print(backgfitcenter,backgheight)
				popt,pcov=cvt(lambda x, A, mu, sigma: gauss2(x,A,mu,sigma,backgheight),xy[int(xlowg):int(xupg)+1],yproj[int(xlowg):int(xupg)+1],p0)

			if int(background) == 1:	
				axp.plot(xy[int(xlow):int(xup)+1],yprojbs[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjY')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				p0=np.array([yprojbs[int(meanfitcenter)],xy[int(meanfitcenter)],xy[2*(int(meanfitcenter)-int(xlowg))]])
				backgheight=np.average(yprojbs[int(backgfitleft):int(backgfitright)])
				print(backgfitcenter,backgheight)
				popt,pcov=cvt(lambda x, A, mu, sigma: gauss2(x,A,mu,sigma,backgheight),xy[int(xlowg):int(xupg)+1],yprojbs[int(xlowg):int(xupg)+1],p0)
			
			axp.plot(xy[int(xlowg):int(xupg)+1],gauss2(xy[int(xlowg):int(xupg)+1],*popt,backgheight),'r',linewidth=0.5)
			area=(np.sqrt(2*(np.pi)))*popt[0]*abs(popt[2])
			print("mean:",popt[1],"sigma:",abs(popt[2]),"area:",area)


#Double Gauss Fit:
	if event.key == 'd':
		print("Perform Gauss Fit:\n")
		print("Set the background left:\n")
		backgfitleft=input()
		print("Set the background right:\n")
		backgfitright=input()
		print("Select the region of interest:\n")
		print("select lower xlim:\n")
		xlowg=input()
		print("select the upper xlim:\n")
		xupg=input()
		print("input the approximate width for first peak")
		print("enter the height of the first peak")
		height1=input()
		print("enter the height of the second peak")
		height2=input()
		print("enter the center of the first peak")
		center1=input()
		print("enter the center of the second peak")
		center2=input()
		print("enter the width of the first peak")
		width1=input()
		print("enter the width of the second peak")
		width2=input()
		backgfitcenter=int(backgfitleft)+int((int(backgfitright)-int(backgfitleft))/2)
		#meanfitcenter=int(xlowg)+int((int(xupg)-int(xlowg))/2)
		axp.clear()
		if int(proj) == 1:
			if int(background) == 0:
				axp.plot(xx[int(xlow):int(xup)+1],xproj[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjX')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection X')
				#Background Height Average:
				backgheight2=np.average(xproj[int(backgfitleft):int(backgfitright)])	
				p0=np.array([int(height1),int(center1),int(width1),int(height2),int(center2),int(width2)])
				print(backgfitcenter,backgheight2)
				popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2: gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight2),xx[int(xlowg):int(xupg)+1],xproj[int(xlowg):int(xupg)+1],p0)

			if int(background) == 1:	
				axp.plot(xx[int(xlow):int(xup)+1],xprojbs[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjX')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				backgheight2=np.average(xprojbs[int(backgfitleft):int(backgfitright)])	
				p0=np.array([int(height1),int(center1),int(width1),int(height2),int(center2),int(width2)])
				print(backgfitcenter,backgheight2)
				popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2: gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight2),xx[int(xlowg):int(xupg)+1],xprojbs[int(xlowg):int(xupg)+1],p0)

#Gauss Fit Results:
			axp.plot(xx[int(xlowg):int(xupg)+1],gaussdoub(xx[int(xlowg):int(xupg)+1],*popt,backgheight2),'r',linewidth=0.5)
			axp.plot(xx[int(xlowg):int(xupg)+1],gausspure(xx[int(xlowg):int(xupg)+1],popt[0],popt[1],popt[2])+backgheight2,'g',linewidth=0.5)
			axp.plot(xx[int(xlowg):int(xupg)+1],gausspure(xx[int(xlowg):int(xupg)+1],popt[3],popt[4],popt[5])+backgheight2,'k',linewidth=0.5)
			area1=(np.sqrt(2*(np.pi)))*popt[0]*abs(popt[2])
			area2=(np.sqrt(2*(np.pi)))*popt[3]*abs(popt[5])
			print("mean1:",popt[1],"sigma1:",abs(popt[2]),"area1:",area1,"mean2:",popt[4],"sigma2:",abs(popt[5]),"area2:",area2)

		if int(proj) == 2:
			if int(background) == 0:
				axp.plot(xy[int(xlow):int(xup)+1],yproj[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjY')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				axp.set_title('Full Projection Y')
				backgheight2=np.average(yproj[int(backgfitleft):int(backgfitright)])	
				p0=np.array([int(height1),int(center1),int(width1),int(height2),int(center2),int(width2)])
				print(backgfitcenter,backgheight2)
				popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2: gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight2),xy[int(xlowg):int(xupg)+1],yproj[int(xlowg):int(xupg)+1],p0)

			if int(background) == 1:	
				axp.plot(xy[int(xlow):int(xup)+1],yprojbs[int(xlow):int(xup)+1],linewidth=0.75,ls='steps',label='ProjY')
				axp.set_title('Gate on '+str(gate)+' keV')
				axp.legend()
				axp.xaxis.set_minor_locator(tck.AutoMinorLocator())
				backgheight2=np.average(yprojbs[int(backgfitleft):int(backgfitright)])	
				p0=np.array([int(height1),int(center1),int(width1),int(height2),int(center2),int(width2)])
				print(backgfitcenter,backgheight2)
				popt,pcov=cvt(lambda x, A1, mu1, sigma1, A2, mu2, sigma2: gaussdoub(x,A1,mu1,sigma1,A2,mu2,sigma2,backgheight2),xy[int(xlowg):int(xupg)+1],yprojbs[int(xlowg):int(xupg)+1],p0)
			

			axp.plot(xy[int(xlowg):int(xupg)+1],gaussdoub(xy[int(xlowg):int(xupg)+1],*popt,backgheight2),'r',linewidth=0.5)
			axp.plot(xy[int(xlowg):int(xupg)+1],gausspure(xy[int(xlowg):int(xupg)+1],popt[0],popt[1],popt[2])+backgheight2,'g',linewidth=0.5)
			axp.plot(xy[int(xlowg):int(xupg)+1],gausspure(xy[int(xlowg):int(xupg)+1],popt[3],popt[4],popt[5])+backgheight2,'k',linewidth=0.5)
			area1=(np.sqrt(2*(np.pi)))*popt[0]*abs(popt[2])
			area2=(np.sqrt(2*(np.pi)))*popt[3]*abs(popt[5])
			print("mean1:",popt[1],"sigma1:",abs(popt[2]),"area1:",area1,"mean2:",popt[4],"sigma2:",abs(popt[5]),"area2:",area2)



	'''
	if event.key == 'n':
		plt.figure(2)
		print("Enter the region to expand\n")
		print("Enter the lower region\n")
		xlowp=input()
		print("Enter the higher region\n")
		xupp=input()
		axp.set_xlim(int(xlowp),int(xupp))
		#axp.autoscale()
		plt.draw()
	'''


#plt.ion()
fig.canvas.mpl_connect('button_press_event',onclick) 
fig.canvas.mpl_connect('key_press_event',onpress)
#figp.canvas.mpl_connect('button_press_event',onclick)
figp.canvas.mpl_connect('key_press_event',onpress)
axp.xaxis.set_minor_locator(tck.AutoMinorLocator())

plt.show()
