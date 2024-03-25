#!/usr/bin/env python3

from struct import *

def DSPsetviewer(fpr,*,offset = 0):
	'''
	Pixie16 DSP viewer
	C. Wibisono
	03/25 '24
	Usage:
	to see Pixie16 DSP parameter from .set file
	Function Argument:
	fpr: file pointer object
	offset: (int) offset value allowed value from 0 up to 23
	'''
	p=Struct("@I")
	fpr.seek(offset*1280*4,0)
	mod1,=p.unpack(fpr.read(4))#mod num
	fpr.seek(-4,1)

	fpr.seek((48*4),1)
	mod1crID,=p.unpack(fpr.read(4))
	mod1slID,=p.unpack(fpr.read(4))
	mod1modID,=p.unpack(fpr.read(4))
	fpr.seek(-12,1)
	fpr.seek(-(48*4),1)
	
	fpr.seek(144*4,1)
	enrise = []
	enflat = []
	trigrise = []
	trigflat = []
	peaksamp = []
	peaksep = []
	CFDtrsh = []
	trigtrsh = []
	for i in range(16):
		temp,=p.unpack(fpr.read(4))
		enrise.append(temp) #slow length
	for j in range(16):
		temp2,=p.unpack(fpr.read(4)) #slow gap
		enflat.append(temp2)
	for k in range(16):
		temp3,=p.unpack(fpr.read(4)) #fast length
		trigrise.append(temp3)
	for l in range(16):
		temp4,=p.unpack(fpr.read(4)) #fast gap
		trigflat.append(temp4)
	for m in range(16):
		temp5,=p.unpack(fpr.read(4))
		peaksamp.append(temp5)
	for n in range(16):
		temp6,=p.unpack(fpr.read(4))
		peaksep.append(temp6)
	for o in range(16):
		temp7,=p.unpack(fpr.read(4))
		CFDtrsh.append(temp7)
	for i in range(16):
		temp8,=p.unpack(fpr.read(4)) #fast treshold
		trigtrsh.append(temp8)
	fpr.seek(-(16*4*8),1)
	fpr.seek(-(144*4),1)
	fpr.seek(352*4,1)
	trlen = []
	for q in range(16):
		temp9,=p.unpack(fpr.read(4))
		trlen.append(temp9)
	fpr.seek(-16*4,1)
	fpr.seek(-(352*4),1)
	fpr.seek(704*4,1)
	qsum0=[]
	qsum1=[]
	qsum2=[]
	qsum3=[]
	qsum4=[]
	qsum5=[]
	qsum6=[]
	qsum7=[]
	for r1 in range(16):
		tempa,=p.unpack(fpr.read(4))
		qsum0.append(tempa)
	 
	for r2 in range(16):
		tempb,=p.unpack(fpr.read(4))
		qsum1.append(tempb)
	
	for r3 in range(16):
		tempc,=p.unpack(fpr.read(4))
		qsum2.append(tempc)
	
	for r4 in range(16):
		tempd,=p.unpack(fpr.read(4))
		qsum3.append(tempd)

	for r5 in range(16):
		tempe,=p.unpack(fpr.read(4))
		qsum4.append(tempe)

	for r6 in range(16):
		tempf,=p.unpack(fpr.read(4))
		qsum5.append(tempf)

	for r7 in range(16):
		tempg,=p.unpack(fpr.read(4))
		qsum6.append(tempg)

	for r8 in range(16):
		temph,=p.unpack(fpr.read(4))
		qsum7.append(temph)
	
	fpr.seek(-(16*8*4),1)
	fpr.seek(-(704*4),1)
	fpr.seek(560*4,1)
	preamptau=[]
	for s in range(16):
		temp10,=p.unpack(fpr.read(4))
		preamptau.append(temp10)
	fpr.seek(-(16*4),1)
	fpr.seek(-(560*4),1)
	paflen = []
	trigdelay = []
	fpr.seek(288*4,1)
	for t1 in range(16):
		temp11,=p.unpack(fpr.read(4))
		paflen.append(temp11)
	for t2 in range(16):
		temp12,=p.unpack(fpr.read(4))
		trigdelay.append(temp12)
	fpr.seek(-(16*2)*4,1)
	fpr.seek(-(288*4),1)

	
	#Output:
	fpr.seek(834*4,1)
	temptimea,= p.unpack(fpr.read(4))
	runtimea = temptimea
	temptimeb, = p.unpack(fpr.read(4)) 
	runtimeb = temptimeb
	fpr.read(4)
	tempcha, = p.unpack(fpr.read(4))
	numevta = tempcha
	tempchb, = p.unpack(fpr.read(4))
	numevtb = tempchb
	fpr.seek(-20,1)
	fpr.seek(-(834*4),1)
	fpr.seek((895*4),1)
	livetimea=[]
	livetimeb=[]
	fastpeaka=[]
	fastpeakb=[]
	for i in range(16):
		templa,=p.unpack(fpr.read(4))
		livetimea.append(templa)
	for j in range(16):
		templb,=p.unpack(fpr.read(4))
		livetimeb.append(templb)
	for k in range(16):
		tempfpa,=p.unpack(fpr.read(4))
		fastpeaka.append(tempfpa)
	for l in range(16):
		tempfpb,=p.unpack(fpr.read(4))
		fastpeakb.append(tempfpb)

	fpr.seek(-(16*4*4),1)
	fpr.seek(-(895*4),1)
	fpr.seek((1055*4),1)
	chnevta=[]
	chnevtb=[]
	for i in range(16):
		tempcha,=p.unpack(fpr.read(4))
		chnevta.append(tempcha)
	for j in range(16):
		tempchb,=p.unpack(fpr.read(4))
		chnevtb.append(tempchb)
	fpr.seek(-(16*2*4),1)
	fpr.seek(-(1055*4),1)
	
	print("mod1:",mod1,"crID:",mod1crID,"slID:",mod1slID,"modID:",mod1modID)
	print("----------------------------------\n")
	
	for i in range(16):
		print("slowrise_t:",enrise[i],"slowflat_t:",enflat[i],"fastrise_t:",trigrise[i],"fastflat_t:",trigflat[i],"fasttrsh:",trigtrsh[i],"cfdtrsh:",CFDtrsh[i],"peak_samp:",peaksamp[i],"peak_sep:",peaksep[i],"trlen:",trlen[i],"paflen:",paflen[i],"trigdelay:",trigdelay[i])
	print("\n")
	print("QDCsumlength:\n")
	for i in range(16):
		print("qsum1:",qsum0[i],"qsum2:",qsum1[i],"qsum3:",qsum2[i],"qsum4:",qsum3[i],"qsum5:",qsum4[i],"qsum6:",qsum5[i],"qsum7:",qsum6[i],"qsum8:",qsum7[i])
	'''
	print("\n")
	print("RunStatistics_______________\n")
	print("runtimea:",runtimea,"runtimeb:",runtimeb,"numevta:",numevta,"numevtb:",numevtb)
	for i in range(16):
		print("chn:",i,"icra:",(fastpeaka[i]/livetimea[i])*100,"icrb:",(fastpeakb[i]/livetimeb[i])*100,"ocra:",(chnevta[i]/runtimea)*100,"ocrb:",(chnevtb[i]/runtimeb)*100)
	'''	


if __name__ == "__main__":
	import sys
	setfile = sys.argv[1]
	with open(setfile, mode ='rb') as f:
		for i in range(24):
			DSPsetviewer(f,offset = i)
