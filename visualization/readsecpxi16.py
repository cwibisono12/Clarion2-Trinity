#!/usr/bin/env python3
#Online Histogram Reader for Ge detector
#C. Wibisono
#05/29/2023


#How to Run:
#./readsecpxi16.py .secfile .calfile dimy dimx xlow xup ylow yup
#Example:
#./readsecpxi16.py Co60.sec /home/clarion/Clarion2-Trinity/param/cal_fsu.ca3 144 8192 100 3000 0 5000
#Above command would display Ge spectra from E=100 to E=3000 with y range starts from 0 to 5000


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import sys 
import os
from scipy.optimize import curve_fit as cvt
#sys.path.insert(1,'/home/clarion/Clarion2-Trinity/visualization/custome')
sys.path.append('custome')
from custome import projmod as p
from custome import fitgui2 as f

filesec=open(sys.argv[1])
filecal=open(sys.argv[2])

dimy=int(sys.argv[3])
dimx=int(sys.argv[4])
xlow=int(sys.argv[5])
xup=int(sys.argv[6])
ylow=int(sys.argv[7])
yup=int(sys.argv[8])

linecal=filecal.readlines()
plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times New Roman']+plt.rcParams['font.serif']

ge=p.clarion(filesec,dimy,dimx,0,dimy-1,0,dimx-1)
ytrp=ge.readparse()

ecal=np.ndarray((dimy,2),dtype=np.float32)
iddet=np.ndarray((dimy,1),dtype=np.int32)

for i in range(0,dimy,1):
  iddet[i]=0
  for j in range(0,2,1):
    ecal[i][j]=0
    


row=0
for line in linecal:
  liner=line.split()
  if liner[0] == '1.0':
    continue
  else:
    iddet[row]=int(liner[0])
    ecal[row][0]=np.float32(liner[1])
    ecal[row][1]=np.float32(liner[2])
    #print('det:',iddet,'intercept:',ecal[row][0],'slope:',ecal[row][1])
    row=row+1
    if row == dimy:
      break
    
A11=ge.project(1,0,0,0,0,0,ytrp)
A12=ge.project(1,0,1,1,0,0,ytrp)
A13=ge.project(1,0,2,2,0,0,ytrp)
A14=ge.project(1,0,3,3,0,0,ytrp)
#A1=gead.project(1,0,1,1,0,0,ytrp)

B21=ge.project(1,0,5,5,0,0,ytrp)
B22=ge.project(1,0,6,6,0,0,ytrp)
B23=ge.project(1,0,7,7,0,0,ytrp)
B24=ge.project(1,0,8,8,0,0,ytrp)
#B2=gead.project(1,0,2,2,0,0,ytrp)

C31=ge.project(1,0,10,10,0,0,ytrp)
C32=ge.project(1,0,11,11,0,0,ytrp)
C33=ge.project(1,0,12,12,0,0,ytrp)
C34=ge.project(1,0,13,13,0,0,ytrp)
#C3=gead.project(1,0,3,3,0,0,ytrp)

D41=ge.project(1,0,16,16,0,0,ytrp)
D42=ge.project(1,0,17,17,0,0,ytrp)
D43=ge.project(1,0,18,18,0,0,ytrp)
D44=ge.project(1,0,19,19,0,0,ytrp)
#D4=gead.project(1,0,4,4,0,0,ytrp)

E51=ge.project(1,0,21,21,0,0,ytrp)
E52=ge.project(1,0,22,22,0,0,ytrp)
E53=ge.project(1,0,23,23,0,0,ytrp)
E54=ge.project(1,0,24,24,0,0,ytrp)
#E5=gead.project(1,0,5,5,0,0,ytrp)

G71=ge.project(1,0,32,32,0,0,ytrp)
G72=ge.project(1,0,33,33,0,0,ytrp)
G73=ge.project(1,0,34,34,0,0,ytrp)
G74=ge.project(1,0,35,35,0,0,ytrp)
#G7=gead.project(1,0,7,7,0,0,ytrp)

H81=ge.project(1,0,37,37,0,0,ytrp)
H82=ge.project(1,0,38,38,0,0,ytrp)
H83=ge.project(1,0,39,39,0,0,ytrp)
H84=ge.project(1,0,40,40,0,0,ytrp)
#H8=gead.project(1,0,8,8,0,0,ytrp)

I91=ge.project(1,0,42,42,0,0,ytrp)
I92=ge.project(1,0,43,43,0,0,ytrp)
I93=ge.project(1,0,44,44,0,0,ytrp)
I94=ge.project(1,0,45,45,0,0,ytrp)
#I9=gead.project(1,0,9,9,0,0,ytrp)

J10_1=ge.project(1,0,48,48,0,0,ytrp)
J10_2=ge.project(1,0,49,49,0,0,ytrp)
J10_3=ge.project(1,0,50,50,0,0,ytrp)
J10_4=ge.project(1,0,51,51,0,0,ytrp)
#J10=gead.project(1,0,10,10,0,0,ytrp)

K11_1=ge.project(1,0,53,53,0,0,ytrp)
K11_2=ge.project(1,0,54,54,0,0,ytrp)
K11_3=ge.project(1,0,55,55,0,0,ytrp)
K11_4=ge.project(1,0,56,56,0,0,ytrp)
#K11=gead.project(1,0,11,11,0,0,ytrp)


x=np.arange(0,8191,1)
x11=np.zeros(8192)
x12=np.zeros(8192)
x13=np.zeros(8192)
x14=np.zeros(8192)
x21=np.zeros(8192)
x22=np.zeros(8192)
x23=np.zeros(8192)
x24=np.zeros(8192)
x31=np.zeros(8192)
x32=np.zeros(8192)
x33=np.zeros(8192)
x34=np.zeros(8192)
x41=np.zeros(8192)
x42=np.zeros(8192)
x43=np.zeros(8192)
x44=np.zeros(8192)
x51=np.zeros(8192)
x52=np.zeros(8192)
x53=np.zeros(8192)
x54=np.zeros(8192)
x71=np.zeros(8192)
x72=np.zeros(8192)
x73=np.zeros(8192)
x74=np.zeros(8192)
x81=np.zeros(8192)
x82=np.zeros(8192)
x83=np.zeros(8192)
x84=np.zeros(8192)
x91=np.zeros(8192)
x92=np.zeros(8192)
x93=np.zeros(8192)
x94=np.zeros(8192)
x101=np.zeros(8192)
x102=np.zeros(8192)
x103=np.zeros(8192)
x104=np.zeros(8192)
x111=np.zeros(8192)
x112=np.zeros(8192)
x113=np.zeros(8192)
x114=np.zeros(8192)

for m in range(0,8192,1):
  if m >= 10 and m<=8190:
    x11[m]=int(ecal[0][1]*x[m]+ecal[0][0])
    x12[m]=int(ecal[1][1]*x[m]+ecal[1][0])
    x13[m]=int(ecal[2][1]*x[m]+ecal[2][0])
    x14[m]=int(ecal[3][1]*x[m]+ecal[3][0])
    x21[m]=int(ecal[5][1]*x[m]+ecal[5][0])
    x22[m]=int(ecal[6][1]*x[m]+ecal[6][0])
    x23[m]=int(ecal[7][1]*x[m]+ecal[7][0])
    x24[m]=int(ecal[8][1]*x[m]+ecal[8][0])
    x31[m]=int(ecal[10][1]*x[m]+ecal[10][0])
    x32[m]=int(ecal[11][1]*x[m]+ecal[11][0])
    x33[m]=int(ecal[12][1]*x[m]+ecal[12][0])
    x34[m]=int(ecal[13][1]*x[m]+ecal[13][0])
    x41[m]=int(ecal[16][1]*x[m]+ecal[16][0])
    x42[m]=int(ecal[17][1]*x[m]+ecal[17][0])
    x43[m]=int(ecal[18][1]*x[m]+ecal[18][0])
    x44[m]=int(ecal[19][1]*x[m]+ecal[19][0])
    x51[m]=int(ecal[21][1]*x[m]+ecal[21][0])
    x52[m]=int(ecal[22][1]*x[m]+ecal[22][0])
    x53[m]=int(ecal[23][1]*x[m]+ecal[23][0])
    x54[m]=int(ecal[24][1]*x[m]+ecal[24][0])
    x71[m]=int(ecal[32][1]*x[m]+ecal[32][0])
    x72[m]=int(ecal[33][1]*x[m]+ecal[33][0])
    x73[m]=int(ecal[34][1]*x[m]+ecal[34][0])
    x74[m]=int(ecal[35][1]*x[m]+ecal[36][0])
    x81[m]=int(ecal[37][1]*x[m]+ecal[37][0])
    x82[m]=int(ecal[38][1]*x[m]+ecal[38][0])
    x83[m]=int(ecal[39][1]*x[m]+ecal[39][0])
    x84[m]=int(ecal[40][1]*x[m]+ecal[40][0])
    x91[m]=int(ecal[42][1]*x[m]+ecal[42][0])
    x92[m]=int(ecal[43][1]*x[m]+ecal[43][0])
    x93[m]=int(ecal[44][1]*x[m]+ecal[44][0])
    x94[m]=int(ecal[45][1]*x[m]+ecal[45][0])
    x101[m]=int(ecal[48][1]*x[m]+ecal[48][0])
    x102[m]=int(ecal[49][1]*x[m]+ecal[49][0])
    x103[m]=int(ecal[50][1]*x[m]+ecal[50][0])
    x104[m]=int(ecal[51][1]*x[m]+ecal[51][0])
    x111[m]=int(ecal[53][1]*x[m]+ecal[53][0])
    x112[m]=int(ecal[54][1]*x[m]+ecal[54][0])
    x113[m]=int(ecal[55][1]*x[m]+ecal[55][0])
    x114[m]=int(ecal[56][1]*x[m]+ecal[56][0])
  else:
    x11[m]=0
    x12[m]=0
    x13[m]=0
    x21[m]=0
    x22[m]=0
    x23[m]=0
    x24[m]=0
    x31[m]=0
    x32[m]=0
    x33[m]=0
    x34[m]=0
    x41[m]=0
    x42[m]=0
    x43[m]=0
    x44[m]=0
    x51[m]=0
    x52[m]=0
    x53[m]=0
    x54[m]=0
    x71[m]=0
    x72[m]=0
    x73[m]=0
    x74[m]=0
    x81[m]=0
    x82[m]=0
    x83[m]=0
    x84[m]=0
    x91[m]=0
    x92[m]=0
    x93[m]=0
    x94[m]=0
    x101[m]=0
    x102[m]=0
    x103[m]=0
    x104[m]=0
    x111[m]=0
    x112[m]=0
    x113[m]=0
    x114[m]=0



fig,ax=plt.subplots(4,1)
fig2,ax2=plt.subplots(4,1)
fig3,ax3=plt.subplots(4,1)
fig4,ax4=plt.subplots(4,1)
fig5,ax5=plt.subplots(4,1)
fig6,ax6=plt.subplots(4,1)
fig7,ax7=plt.subplots(4,1)
fig8,ax8=plt.subplots(4,1)
fig9,ax9=plt.subplots(4,1)
fig10,ax10=plt.subplots(4,1)

for i in range(0,4,1):
	ax[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax2[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax2[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax2[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax3[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax3[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax3[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax4[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax4[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax4[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax5[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax5[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax5[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax6[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax6[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax6[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax7[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax7[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax7[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax8[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax8[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax8[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax9[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax9[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax9[i].xaxis.set_minor_locator(tck.AutoMinorLocator())
	ax10[i].tick_params(direction='in',axis='both',which='major',bottom='False',left='False',top='False',right='False',length=9,width=0.75)
	ax10[i].tick_params(direction='in',axis='both',which='minor',bottom='False',left='False',top='False',right='False',length=6,width=0.75)
	ax10[i].xaxis.set_minor_locator(tck.AutoMinorLocator())

line11,=ax[0].plot(x11,A11,linewidth=0.85,ls='steps-post',label='A1-Black')
line12,=ax[1].plot(x12,A12,linewidth=0.85,ls='steps-post',label='A1-Blue')
line13,=ax[2].plot(x13,A13,linewidth=0.85,ls='steps-post',label='A1-Gree')
line14,=ax[3].plot(x14,A14,linewidth=0.85,ls='steps-post',label='A1-Red')
#ax[4].plot(x[xlow:xup],A1[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='A1-addback')
line101,=ax10[0].plot(x111,K11_1,linewidth=0.85,ls='steps-post',label='K11-Black')
line102,=ax10[1].plot(x112,K11_2,linewidth=0.85,ls='steps-post',label='K11-Blue')
line103,=ax10[2].plot(x113,K11_3,linewidth=0.85,ls='steps-post',label='K11-Green')
line104,=ax10[3].plot(x114,K11_4,linewidth=0.85,ls='steps-post',label='K11-Red')
#ax10[4].plot(x[xlow:xup],K11[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='K11-addback')
line21,=ax2[0].plot(x21,B21,linewidth=0.85,ls='steps-post',label='B2-Black')
line22,=ax2[1].plot(x22,B22,linewidth=0.85,ls='steps-post',label='B2-Blue')
line23,=ax2[2].plot(x23,B23,linewidth=0.85,ls='steps-post',label='B2-Green')
line24,=ax2[3].plot(x24,B24,linewidth=0.85,ls='steps-post',label='B2-Red')
#ax2[4].plot(x[xlow:xup],B2[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='B2-addback')
line31,=ax3[0].plot(x31,C31,linewidth=0.85,ls='steps-post',label='C3-Black')
line32,=ax3[1].plot(x32,C32,linewidth=0.85,ls='steps-post',label='C3-Blue')
line33,=ax3[2].plot(x33,C33,linewidth=0.85,ls='steps-post',label='C3-Gree')
line34,=ax3[3].plot(x34,C34,linewidth=0.85,ls='steps-post',label='C3-Red')
#ax3[4].plot(x[xlow:xup],C3[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='C3-addback')
line41,=ax4[0].plot(x41,D41,linewidth=0.85,ls='steps-post',label='D4-Black')
line42,=ax4[1].plot(x42,D42,linewidth=0.85,ls='steps-post',label='D4-Blue')
line43,=ax4[2].plot(x43,D43,linewidth=0.85,ls='steps-post',label='D4-Green')
line44,=ax4[3].plot(x44,D44,linewidth=0.85,ls='steps-post',label='D4-Red')
#ax4[4].plot(x[xlow:xup],D4[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='D4-addback')
line51,=ax5[0].plot(x51,E51,linewidth=0.85,ls='steps-post',label='E5-Black')
line52,=ax5[1].plot(x52,E52,linewidth=0.85,ls='steps-post',label='E5-Blue')
line53,=ax5[2].plot(x53,E53,linewidth=0.85,ls='steps-post',label='E5-Green')
line54,=ax5[3].plot(x54,E54,linewidth=0.85,ls='steps-post',label='E5-Red')
#ax5[4].plot(x[xlow:xup],E5[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='E5-addback')
line61,=ax6[0].plot(x71,G71,linewidth=0.85,ls='steps-post',label='G7-Black')
line62,=ax6[1].plot(x72,G72,linewidth=0.85,ls='steps-post',label='G7-Blue')
line63,=ax6[2].plot(x73,G73,linewidth=0.85,ls='steps-post',label='G7-Green')
#ax6[4].plot(x[xlow:xup],G7[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='G7-addback')
line64,=ax6[3].plot(x74,G74,linewidth=0.85,ls='steps-post',label='G7-Red')
line71,=ax7[0].plot(x81,H81,linewidth=0.85,ls='steps-post',label='H8-Black')
line72,=ax7[1].plot(x82,H82,linewidth=0.85,ls='steps-post',label='H8-Blue')
line73,=ax7[2].plot(x83,H83,linewidth=0.85,ls='steps-post',label='H8-Green')
line74,=ax7[3].plot(x84,H84,linewidth=0.85,ls='steps-post',label='H8-Red')
#ax7[4].plot(x[xlow:xup],H8[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='H8-addback')
line81,=ax8[0].plot(x91,I91,linewidth=0.85,ls='steps-post',label='I9-Black')
line82,=ax8[1].plot(x92,I92,linewidth=0.85,ls='steps-post',label='I9-Blue')
line84,=ax8[3].plot(x94,I94,linewidth=0.85,ls='steps-post',label='I9-Red')
#ax8[4].plot(x[xlow:xup],I9[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='I9-addback')
line83,=ax8[2].plot(x93,I93,linewidth=0.85,ls='steps-post',label='I9-Green')
line91,=ax9[0].plot(x101,J10_1,linewidth=0.85,ls='steps-post',label='J10-Black')
line92,=ax9[1].plot(x102,J10_2,linewidth=0.85,ls='steps-post',label='J10-Blue')
line93,=ax9[2].plot(x103,J10_3,linewidth=0.85,ls='steps-post',label='J10-Green')
line94,=ax9[3].plot(x104,J10_4,linewidth=0.85,ls='steps-post',label='J10-Red')
#axi9[4].plot(x[xlow:xup],J10[xlow:xup],linewidth=0.85,ls='steps-post',color='red',label='J10-addback')


for i in range(0,4,1):
	ax[i].legend()
	ax[i].set_xlim(left=xlow,right=xup,auto=False)
	ax[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax2[i].legend()
	ax2[i].set_xlim(left=xlow,right=xup,auto=False)
	ax2[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax3[i].legend()
	ax3[i].set_xlim(left=xlow,right=xup,auto=False)
	ax3[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax4[i].legend()
	ax4[i].set_xlim(left=xlow,right=xup,auto=False)
	ax4[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax5[i].legend()
	ax5[i].set_xlim(left=xlow,right=xup,auto=False)
	ax5[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax6[i].legend()
	ax6[i].set_xlim(left=xlow,right=xup,auto=False)
	ax6[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax7[i].legend()
	ax7[i].set_xlim(left=xlow,right=xup,auto=False)
	ax7[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax8[i].legend()
	ax8[i].set_xlim(left=xlow,right=xup,auto=False)
	ax8[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax9[i].legend()
	ax9[i].set_xlim(left=xlow,right=xup,auto=False)
	ax9[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax10[i].legend()
	ax10[i].set_xlim(left=xlow,right=xup,auto=False)
	ax10[i].set_ylim(bottom=ylow,top=yup,auto=False)
	ax[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax2[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax3[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax4[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax5[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax6[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax7[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax8[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax9[i].set_ylabel('counts/keV',style='normal',fontweight='bold')
	ax10[i].set_ylabel('counts/keV',style='normal',fontweight='bold')


ax[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax2[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax3[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax4[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax5[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax6[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax7[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax8[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax9[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
ax10[3].set_xlabel('Energy (keV)',style='normal',fontweight='bold')
fig.suptitle("A1")
fig2.suptitle("B2")
fig3.suptitle("C3")
fig4.suptitle("D4")
fig5.suptitle("E5")
fig6.suptitle("G7")
fig7.suptitle("H8")
fig8.suptitle("I9")
fig9.suptitle("J10")
fig10.suptitle("K11")

p11=f.fitgui(line11)
p12=f.fitgui(line12)
p13=f.fitgui(line13)
p14=f.fitgui(line14)

p21=f.fitgui(line21)
p22=f.fitgui(line22)
p23=f.fitgui(line23)
p24=f.fitgui(line24)

p31=f.fitgui(line31)
p32=f.fitgui(line32)
p33=f.fitgui(line33)
p34=f.fitgui(line34)

p41=f.fitgui(line41)
p42=f.fitgui(line42)
p43=f.fitgui(line43)
p44=f.fitgui(line44)

p51=f.fitgui(line51)
p52=f.fitgui(line52)
p53=f.fitgui(line53)
p54=f.fitgui(line54)

p61=f.fitgui(line61)
p62=f.fitgui(line62)
p63=f.fitgui(line63)
p64=f.fitgui(line64)


p71=f.fitgui(line71)
p72=f.fitgui(line72)
p73=f.fitgui(line73)
p74=f.fitgui(line74)

p81=f.fitgui(line81)
p82=f.fitgui(line82)
p83=f.fitgui(line83)
p84=f.fitgui(line84)

p91=f.fitgui(line91)
p92=f.fitgui(line92)
p93=f.fitgui(line93)
p94=f.fitgui(line94)

p101=f.fitgui(line101)
p102=f.fitgui(line102)
p103=f.fitgui(line103)
p104=f.fitgui(line104)

p11.connect()
p12.connect()
p13.connect()
p14.connect()

p21.connect()
p22.connect()
p23.connect()
p24.connect()

p31.connect()
p32.connect()
p33.connect()
p34.connect()

p41.connect()
p42.connect()
p43.connect()
p44.connect()

p51.connect()
p52.connect()
p53.connect()
p54.connect()

p61.connect()
p62.connect()
p63.connect()
p64.connect()

p71.connect()
p72.connect()
p73.connect()
p74.connect()

p81.connect()
p82.connect()
p83.connect()
p84.connect()

p91.connect()
p92.connect()
p93.connect()
p94.connect()

p101.connect()
p102.connect()
p103.connect()
p104.connect()



plt.show()
filesec.close()
filecal.close()      
