#!/usr/bin/env python3
from struct import *
from dataclasses import dataclass


@dataclass
class pxi16:
	chn: int
	sln: int
	crn: int
	iddet: int
	hlen: int
	elen: int
	fcode: int
	time: int
	ctime: int
	ctimef: int
	energy: int
	trlen: int
	trwlen: int
	extra: int
	#esum: list[int]
	#qsum: list[int]
	#tr: list[int]

def pxi16file(fpr):
	'''
	Original Pixie16 Data Format
	C. Wibisono
	02/16 '24	
	Usage:
	To parse the .evt file coming from original pixie16 data structure per detector hit. 
	Function Argument:
	fpr: file pointer object
	Return value: (either -1, or pxi16data type (see pxi16 dataclasses above)
	-1: end of file
	pxi16data: pixie16 data type object
	'''
	buff1=fpr.read(4)
	if buff1 == b'':
		return -1 #break from the iteration over channel
	p=Struct("@i")
	buff1int,=p.unpack(buff1)
	chn =int(buff1int & 0xF)
	sln=int((buff1int & 0xF0) >> 4)
	crn=int((buff1int & 0XF00) >> 8)
	hlen=int((buff1int & 0x1F000) >> 12)
	elen=int((buff1int & 0x7FFE0000) >> 17)
	fcode=int((buff1int & 0x80000000) >> 31)
	buff2=fpr.read(4)
	buff3=fpr.read(4)
	buff2int,=p.unpack(buff2)
	buff3int,=p.unpack(buff3)
	time=int(((buff3int & 0xFFFF) << 32) + buff2int)			
	ctime=int((buff3int & 0x7FFF0000) >> 16)
	ctimef=int((buff3int & 0x80000000) >> 31)
	buff4=fpr.read(4)
	buff4int,=p.unpack(buff4)
	tempenergy=int (buff4int & 0xFFFF)
	tempen=float(tempenergy/2.)
	energy=int(tempen)
	trlen=int ((buff4int & 0x7FFF0000) >> 16)
	trwlen=int(trlen/2.)  
	extra=int((buff4int & 0x80000000) >> 31)
	iddet=crn*13*16+(sln-2)*16+chn	
	if hlen == 4 and trwlen == 0:
		ans=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
		return ans  #continue iteration over channel
	#f.seek(-16,1)
	esum=[]
	qsum=[]
	if hlen == 8 or hlen == 16:
		for i in range(0,4,1):
			buff=fpr.read(4)
			temp,=p.unpack(buff)
			esum.append(temp)
		if hlen == 16:
			for j in range(0,8,1):
				buff=fpr.read(4)
				temp,=p.unpack(buff)	
				qsum.append(temp)

	if hlen == 12:
		for i in range(0,8,1):
			buff=fpr.read(4)
			temp,=p.unpack(buff)
			qsum.append(temp)

		ans=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
		return ans

	tr=[]
	for i in range(hlen,elen,1):
		buff=fpr.read(4)
		buffint,=p.unpack(buff)				
		if trwlen !=0:
			temp1=int(buffint & 0x3FFF)
			temp2=int((buffint >> 16) & 0x3FFF)
			tr.append(temp1)
			tr.append(temp2)
			
			
	if hlen == 4 and trwlen !=0:
		temp1=0
		for i in range(0,31,1):
			temp1=temp1+tr[i]
		qsum.append(temp1)	
		temp2=0
		for i in range(31,60,1):
			temp2=temp2+tr[i]
		qsum.append(temp2)
		temp3=0
		for i in range(60,75,1):
			temp3=temp3+tr[i]
		qsum.append(temp3)
		temp4=0
		for i in range(75,95,1):
			temp4=temp4+tr[i]
		qsum.append(temp4)
		temp5=0
		for i in range(95,105,1):
			temp5=temp5+tr[i]
		qsum.append(temp5)
		temp6=0
		for i in range(105,160,1):
			temp6=temp6+tr[i]
		qsum.append(temp6)
		temp7=0
		for i in range(160,175,1):
			temp7=temp7+tr[i]
		qsum.append(temp7)
		temp8=0
		for i in range(175,200,1):
			temp8=temp8+tr[i]
		qsum.append(temp8)

		print("chn:",chn,"sln:",sln,"crn:",crn,"hlen:",hlen,"elen:",elen,"energy:",energy)
		ans=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
		return ans 

def nsclpxi16file(fpr):
	'''
	NSCL/FRIB Pixie16 Data Format
	C. Wibisono
	02/16 '24
	Usage:
	To parse the .evt file coming from NSCL/FRIB DAQ. This function parser would return pxi16data 		object per event.
	Function Argument(s):
	fpr: file pointer object
	Return value: either -1, 0 or pxi16 data type
	-1: end of file
	0: skip for non-physics ring
	pxi16dataobj[sevtmult]: list of pxi16data object. The dimension would denote number of hits fo		r each event.
	'''
	p=Struct("@i")
	buffring=fpr.read(4)
	if buffring == b'':
		return -1
	rihsize,=p.unpack(buffring)
	bufftype=fpr.read(4)
	rihtype,=p.unpack(bufftype)
	if rihtype != 30:
		fpr.seek(rihsize-8,1)
		return 0
	buffsize=fpr.read(4)
	ribhsize,=p.unpack(buffsize)
	buffht=fpr.read(8)
	ribht,=unpack("@q",buffht)
	buffsid=fpr.read(4)
	ribhsid,=p.unpack(buffsid)
	buffhbt=fpr.read(4)
	ribhbt,=p.unpack(buffhbt)
	buffbhsize=fpr.read(4)
	ribhbsize,=p.unpack(buffbhsize)
	temporary=ribhbsize-4
	#iterate for each fragment for a given ring
	sevtmult=0
	pxi16obj=[]
	while(temporary>0):
		fpr.read(8) #frtmp
		fpr.read(4) #sid
		fpr.read(4) #payload
		fpr.read(4) #bt
		fpr.read(4) #rihsize
		fpr.read(4) #rihtype
		fpr.read(4) #ribhsize
		fpr.read(8) #temp
		fpr.read(4) #sid2
		fpr.read(4) #bt2
		buffbsize=fpr.read(4) #bsize
		bsize,=p.unpack(buffbsize)
		buffdevice=fpr.read(4)
		if buffdevice == b'':
			return -1
		temporary=temporary-(48+2*bsize)
		buff1=fpr.read(4)
		if buff1 == b'':
			return -1
		buff1int,=p.unpack(buff1)
		chn =int(buff1int & 0xF)
		sln=int((buff1int & 0xF0) >> 4)
		crn=int((buff1int & 0XF00) >> 8)
		hlen=int((buff1int & 0x1F000) >> 12)
		elen=int((buff1int & 0x7FFE0000) >> 17)
		fcode=int((buff1int & 0x80000000) >> 31)
		buff2=fpr.read(4)
		buff3=fpr.read(4)
		buff2int,=p.unpack(buff2)
		buff3int,=p.unpack(buff3)
		time=int(((buff3int & 0xFFFF) << 32) + buff2int)			
		ctime=int((buff3int & 0x7FFF0000) >> 16)
		ctimef=int((buff3int & 0x80000000) >> 31)
		buff4=fpr.read(4)
		buff4int,=p.unpack(buff4)
		tempenergy=int (buff4int & 0xFFFF)
		tempen=float(tempenergy/2.)
		energy=int(tempen)
		trlen=int ((buff4int & 0x7FFF0000) >> 16)
		trwlen=int(trlen/2.)  
		extra=int((buff4int & 0x80000000) >> 31)
		iddet=crn*13*16+(sln-2)*16+chn	
			
		if hlen == 4 and trwlen == 0:
			sevtmult=sevtmult+1
			#print("chn:",chn,"sln:",sln,"crn:",crn,"hlen:",hlen,"elen:",elen,"energy:",energy)
			temp=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
			pxi16obj.append(temp)
			continue
			#f.seek(-16,1)
		esum=[]
		qsum=[]
		if hlen == 8 or hlen == 16:
			for i in range(0,4,1):
				buff=fpr.read(4)
				temp,=p.unpack(buff)
				esum.append(temp)
			if hlen == 16:
				for j in range(0,8,1):
					buff=fpr.read(4)
					temp,=p.unpack(buff)	
					qsum.append(temp)

		if hlen == 12:
			for i in range(0,8,1):
				buff=fpr.read(4)
				temp,=p.unpack(buff)
				qsum.append(temp)

			
		tr=[]
		for i in range(hlen,elen,1):
			buff=fpr.read(4)
			buffint,=p.unpack(buff)
				
			if trwlen !=0:
				temp1=int(buffint & 0x3FFF)
				temp2=int((buffint >> 16) & 0x3FFF)
				tr.append(temp1)
				tr.append(temp2)
			
				'''
				if hlen == 4 and trwlen !=0:
					temp1=0
					for i in range(0,31,1):
						temp1=temp1+tr[i]
					qsum.append(temp1)	
					temp2=0
					for i in range(31,60,1):
						temp2=temp2+tr[i]
					qsum.append(temp2)
					temp3=0
					for i in range(60,75,1):
						temp3=temp3+tr[i]
					qsum.append(temp3)
					temp4=0
					for i in range(75,95,1):
						temp4=temp4+tr[i]
					qsum.append(temp4)
					temp5=0
					for i in range(95,105,1):
						temp5=temp5+tr[i]
					qsum.append(temp5)
					temp6=0
					for i in range(105,160,1):
						temp6=temp6+tr[i]
					qsum.append(temp6)
					temp7=0
					for i in range(160,175,1):
						temp7=temp7+tr[i]
					qsum.append(temp7)
					temp8=0
					for i in range(175,200,1):
						temp8=temp8+tr[i]
					qsum.append(temp8)
				'''
		sevtmult=sevtmult+1
		#print("chn:",chn,"sln:",sln,"crn:",crn,"hlen:",hlen,"elen:",elen,"energy:",energy)
		temp=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
		pxi16obj.append(temp)
	
	return pxi16obj

if __name__ == "__main__":
	'''
	Below is an example of how to use parsegen and pxi16 parser module
	to create an id-Energy Matrix from the raw data.
	'''
	import sys
	from parsegen import matwrite
	import numpy as np
	filename=sys.argv[1]
	matfilename=sys.argv[2]
	overwrite=int(sys.argv[3]) #1 to overwrite, 0 to append
	
	idenmat=np.ndarray(shape=(416,8192),dtype=np.int32) #idvsenergy 2D mat (int 4 bytes)
	
	#open the raw data
	with open(filename,mode='rb') as f:
		while(1):
			temp=nsclpxi16file(f)
			if temp == -1:
				break
			elif temp == 0:
				continue
			else:
				sevtmult=len(temp) #event multiplicity
				for i in range(0,sevtmult,1):
					print("sevtmult:",sevtmult,"crn:",temp[i].crn,"sln:",temp[i].sln,"chn:",temp[i].chn,temp[i].iddet)
					if temp[i].energy >=0 and temp[i].energy < 8192:
						idenmat[temp[i].iddet][temp[i].energy]=idenmat[temp[i].iddet][temp[i].energy]+1
					
	#Write matrix into a file
	matwrite(matfilename,dimy=416,dimx=8192,arr=idenmat,overwrite=overwrite)
