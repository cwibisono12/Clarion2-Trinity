#!/usr/bin/env python3

#./2Dvis.py file 5000 5000 0 4999 0 4999 
#C. Wibisono



#Program for Displaying 2D Histogram + Drawing Banana Cuts + Making Basic Projections:

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.colors import LogNorm

print("2D visualization with the projection\n")
print("Hit g to make Projection")
print("Hit e to expand axis")
print("Hit f to generate coordinates useful for making banana gates")
print("click mouse to draw banana gates")
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
#axp=[]
#Prepared the array for y and x:
yproj=np.zeros(int(sys.argv[2]),dtype=np.int32) #Project onto Y axis
xproj=np.zeros(int(sys.argv[3]),dtype=np.int32) #Project onto X axis
xx=np.zeros(int(sys.argv[3]),dtype=np.int32)
xy=np.zeros(int(sys.argv[2]),dtype=np.int32)
#x=np.arange(0,int(sys.argv[2]),1)
#x=np.zeros(100)
#y=np.zeros(100)
#for i in range(100):
#	x[i]=i
#	y[i]=2*x[i]

#ax.plot(x,y)

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
		global axp,xx,xy,proj
		if int(waxis) == 1:
			print("Enter the lower xlim\n")
			xlow=input()
			print("Enter the upper xlim\n")
			xup=input()
			print("Enter the lower ylim\n")
			ylow=input()
			print("Enter the upper ylim\n")
			yup=input()
			ax.set_xlim(int(xlow),int(xup))
			ax.set_ylim(int(ylow),int(yup))
		if int(waxis) == 2: 
			print("Enter the lower limit\n")
			xlow=input()
			print("Enter the upper limit\n")
			xup=input()
			axp.clear()
			axp.set_xlim(int(xlow),int(xup))
			if int(proj)==1:	
				axp.plot(xx[int(xlow):int(xup)],xproj[int(xlow):int(xup)],linewidth=0.75,ls='steps',label='ProjX')
				axp.legend()
			if int(proj)==2:
				axp.plot(xy[int(xlow):int(xup)],yproj[int(xlow):int(xup)],linewidth=0.75,ls='steps',label='ProjY')
				axp.legend()
			#axp.autoscale()
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
		print("Enter the region to Project:\n")
		print("Enter the lower region:\n")
		low=input()
		print("Enter the higher region:\n")
		high=input()
		#global axp
		#Projection on X axis: Make a gate on Y axis then project on X axis:
		axp.clear()
		if int(proj) == 1:
			#global x
			xx=np.arange(0,int(sys.argv[3]),1)
			for i in range(0,int(sys.argv[2]),1):
				for j in range(int(low),int(high)+1,1):
					xproj[i]=xproj[i]+ytrp[j,i]	
			axp.plot(xx,xproj,linewidth=0.75,ls='steps',label='ProjX')
			axp.legend()
		#Projection on Y axis: Make a gate on X axis then project on Y axis:
		if int(proj) == 2:
			#global x
			xy=np.arange(0,int(sys.argv[2]),1)
			for i in range(int(low),int(high)+1,1):
				for j in range(0,int(sys.argv[2]),1):
					yproj[i]=yproj[i]+ytrp[i,j]
			axp.plot(xy,yproj,linewidth=0.75,ls='steps',label='ProjY')
			axp.legend()
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
fig.canvas.mpl_connect('button_press_event',onclick) 
fig.canvas.mpl_connect('key_press_event',onpress)
figp.canvas.mpl_connect('button_press_event',onclick)
figp.canvas.mpl_connect('key_press_event',onpress)
plt.show()
