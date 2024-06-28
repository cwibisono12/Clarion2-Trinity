#!/usr/bin/env python3

from struct import *


def pxi16evfilter(fpr, *, timebuild=190):
	'''
	Pixie16 List-Mode Data Reduction Module.
	C. Wibisono
	06/27 '24	
	Usage:
	To filter pixie16 list mode data. Current filter requirement is based on the data that have gamma. Any event that does not have gamma hits will be filtered out. This function can be modified based on specific user's requirement.

	Function Argument(s):
	fpr: file pointer object
	timebuild: time event window in the unit of 100 MHz (10 ns).
	Return value: -1, 0, 1  or pxi16data type (see pxi16 dataclasses above)
	-1: error indicating that the file is not time sorted
	0: end of file
	1: flag indicating that the event does not meet the requirement and will be filtered out.

	pxi16raw[sevtmult]: pixie16 data type object. The number of dimension for pxi16data object denotes the number of detector hits within an event.
	'''
	evtime = -1
	tdiff = -1
	pxi16raw = []
	iddetlist = []
	p=Struct("@I")

	while(1):
		buff1=fpr.read(4)
		if buff1 == b'':
			return 0
		buff1int,=p.unpack(buff1)
		chn =buff1int & 0xF
		sln=(buff1int & 0xF0) >> 4
		crn=(buff1int & 0XF00) >> 8
		hlen=(buff1int & 0x1F000) >> 12
		elen=(buff1int & 0x7FFE0000) >> 17
		fcode=(buff1int & 0x80000000) >> 31
		buff2=fpr.read(4)
		buff3=fpr.read(4)
		buff2int,=p.unpack(buff2)
		buff3int,=p.unpack(buff3)
		time=((buff3int & 0xFFFF) << 32) + (buff2int)	
		ctime=(buff3int & 0x7FFF0000) >> 16
				
		#ctime=(buff3int & 0xFFFF0000) >> 16
		
		ctimef=(buff3int & 0x80000000) >> 31
		buff4=fpr.read(4)
		buff4int,=p.unpack(buff4)
		tempenergy= buff4int & 0xFFFF
		tempen=float(tempenergy/2.)
		energy=int(tempen)
		trlen= (buff4int & 0x7FFF0000) >> 16
		trwlen=int(trlen/2.)  
		extra= (buff4int & 0x80000000) >> 31
		iddet=crn*13*16+(sln-2)*16+chn	
		
	
		if evtime == -1:
			evtime = time
			tdiff = 0
		else:
			tdiff = time - evtime
			if tdiff < 0:
				print("file is not time sorted.\n")
				return -1
		
		#rewind the event if it exceeds the time window
		if tdiff > timebuild:
			fpr.seek(-16,1)
			break
		
		#Record all detector id for each event:
		iddetlist.append(iddet)		
		
		if hlen == 4 and trwlen == 0:
			pxi16raw.append([buff1int, buff2int, buff3int, buff4int])
			continue  #continue iteration over channel
	
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
		trraw=[]
		for i in range(hlen,elen,1):
			buff=fpr.read(4)
			buffint,=p.unpack(buff)				
			trraw.append(buffint)

			if trwlen !=0:
				temp1=buffint & 0x3FFF
				temp2=(buffint >> 16) & 0x3FFF
				tr.append(temp1)
				tr.append(temp2)
			
		

		if trwlen == 0:
			if hlen == 8:
				pxi16raw.append([buff1int,buff2int,buff3int,buff4int,esum])	
			if hlen == 12:
				pxi16raw.append([buff1int,buff2int,buff3int,buff4int,qsum])
			if hlen == 16:
				pxi16raw.append([buff1int,buff2int,buff3int,buff4int,esum,qsum])
		else:
			if hlen == 8:
				pxi16raw.append([buff1int,buff2int,buff3int,buff4int,esum,trraw])
			if hlen == 12:
				pxi16raw.append([buff1int,buff2int,buff3int,buff4int,qsum,trraw])
			if hlen == 16:
				pxi16raw.append([buff1int,buff2int,buff3int,buff4int,esum,qsum,trraw])


	num = len(iddetlist)
	Ge = 0
	GAGG = 0

	for k in range(num):
		if iddetlist[k] >= 0 and iddetlist[k] <= 57:
			Ge = Ge + 1
		if iddetlist[k] >= 64 and iddetlist[k] <=191:
			GAGG = GAGG + 1

	if Ge > 0:
		return pxi16raw
	if Ge == 0:
		return 1

def pxi16evwrite(fout, pxi16raw):
	'''Function to write Pixie16 Data Format into a file.
	C. Wibisono
	06/27 '24
	Parameter(s):
	fout: file pointer object to write into a file.
	pxi16raw: Pixie16 pointer object data format.
	Usage:
	Along with the pxi16evfilter function, user can write the filtered version of pxi16 list mode data while preserving the structure of the Pixie16 List Mode Data.
	'''

	sevtmult = len(pxi16raw)
	p = Struct("@I")
	
	for i in range(sevtmult):
		#Write 4 words header:
		num = len(pxi16raw[i])
		fout.write(p.pack(pxi16raw[i][0]))
		fout.write(p.pack(pxi16raw[i][1]))
		fout.write(p.pack(pxi16raw[i][2]))
		fout.write(p.pack(pxi16raw[i][3]))
		
		#More than 4 words header:
		if num == 5:
			dim=len(pxi16raw[i][4])
			if dim == 4:
				for k in range(4):
					fout.write(p.pack(pxi16raw[i][4][k]))
			if dim == 8:
				for m in range(8):
					fout.write(p.pack(pxi16raw[i][4][m]))

		if num == 6:
			dim4 = len(pxi16raw[i][4])
			dim5 = len(pxi16raw[i][5])
			for k1 in range(dim4):
				fout.write(p.pack(pxi16raw[i][4][k1]))
			for k2 in range(dim5):
				fout.write(p.pack(pxi16raw[i][5][k2]))

		if num == 7:
			dim4 = len(pxi16raw[i][4])
			dim5 = len(pxi16raw[i][5])
			dim6 = len(pxi16raw[i][6])
			for k1 in range(dim4):
				fout.write(p.pack(pxi16raw[i][4][k1]))
			for k2 in range(dim5):
				fout.write(p.pack(pxi16raw[i][5][k2]))
			for k3 in range(dim6):
				fout.write(p.pack(pxi16raw[i][6][k3]))
	



if __name__ == "__main__":
	'''
	Below is an example of how to use pxi16evfilter and pxi16evwrite module to filter pxi16 list mode data.
	'''
	import sys
	
	filein=sys.argv[1] #Original Pixie16 Time-Sorted List-Mode File
	fileout=sys.argv[2] #Filtered Pixie16 List-Mode File
	
	
	with open(fileout, mode='wb') as fout:
		#Open the raw data:
		with open(filein,mode='rb') as fin:
			while(1):
				temp=pxi16evfilter(fin)
				if temp == -1 or temp == 0:
					break
				elif temp == 1:
					continue
				else:
					pxi16evwrite(fout, temp)
