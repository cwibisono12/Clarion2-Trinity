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


def pxi16evread(fpr,*,timebuild=190):
	'''
	Original Pixie16 Data Format
	C. Wibisono
	02/16 '24	
	Usage:
	To parse the .evt file coming from original pixie16 data structure and build one-event. File needs to be time sorted before using this function.
	Function Argument(s):
	fpr: file pointer object
	timebuild: time event window in the unit of 100 MHz (10 ns).
	Return value: -1, 0  or pxi16data type (see pxi16 dataclasses above)
	-1: error indicating that the file is not time sorted
	0: end of file
	pxi16data[sevtmult]: pixie16 data type object. The number of dimension for pxi16data object denotes the number of detector hits within an event.
	'''
	evtime=-1
	tdiff=-1
	pxi16obj=[]
	p=Struct("@I")
	while(1):
		buff1=fpr.read(4)
		if buff1 == b'':
			return 0
			#break #break from the iteration over channel
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

		if hlen == 4 and trwlen == 0:
			temp=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
			pxi16obj.append(temp)
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
		for i in range(hlen,elen,1):
			buff=fpr.read(4)
			buffint,=p.unpack(buff)				
		
			if trwlen !=0:
				temp1=buffint & 0x3FFF
				temp2=(buffint >> 16) & 0x3FFF
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

		
		temp=pxi16(chn,sln,crn,iddet,hlen,elen,fcode,time,ctime,ctimef,energy,trlen,trwlen,extra)
		pxi16obj.append(temp)

	return pxi16obj


def ev5read(f):
        '''
        .ev5 file format
        C. Wibisono
        02/14 '24
        Usage:
        To parse the .ev5 file.
	Parameter(s):
	f: file pointer object
	Return:
	List [ID, Ge (dict)]
        '''
        p=Struct("@b")
        q=Struct("@h")
        buff1=f.read(1)
        if buff1 == b'':
                return -1
        gaggvalidp,=p.unpack(buff1)
        buff2=f.read(1)
        if buff2 == b'':
                return -1
        gaggvalida,=p.unpack(buff2)
        buff3=f.read(1)
        if buff3 == b'':
                return -1
        gmult,=p.unpack(buff3)
	
        if gaggvalidp != 0 or gaggvalida !=0:
                buff4=f.read(2)
                if buff4 == b'':
                        return -1
                    
                Ex,=q.unpack(buff4)
                ID = [gaggvalidp, gaggvalida, gmult, Ex]

        else:
                ID = [gaggvalidp, gaggvalida, gmult]
		

        Ge = {}
        for i in range(0,gmult,1):
                xid=[]
                xevalid=[]
                gid,=q.unpack(f.read(2)) #Clover ID
                gaold,=q.unpack(f.read(2)) #Old Angle
                ganew,=q.unpack(f.read(2)) #New Angle
                gae,=q.unpack(f.read(2)) #Energy 
                gaetrue,=q.unpack(f.read(2)) #Edopp
                gavalidp,=p.unpack(f.read(1)) #validp
                gavalidnp,=p.unpack(f.read(1)) #validnp
                gatime,=q.unpack(f.read(2)) #xmult
                gaxmult,=p.unpack(f.read(1))
                Ge[str(gid)] = [gaold, ganew, gae, gaetrue, gavalidp, gavalidnp, gatime, gaxmult]
                temp1a,=p.unpack(f.read(1)) #xid
                temp1b,=p.unpack(f.read(1))
                temp1c,=p.unpack(f.read(1))
                temp1d,=p.unpack(f.read(1))
                f.read(1)
                temp2a,=q.unpack(f.read(2)) #xevalid
                temp2b,=q.unpack(f.read(2))
                temp2c,=q.unpack(f.read(2))
                temp2d,=q.unpack(f.read(2))
                xid.append(temp1a)
                xid.append(temp1b)
                xid.append(temp1c)
                xid.append(temp1d)
                xevalid.append(temp2a)
                xevalid.append(temp2b)
                xevalid.append(temp2c)
                xevalid.append(temp2d)
                Ge[str(gid)].append(xid)
                Ge[str(gid)].append(xevalid)
                #print("p: ",int(gaggvalidp),"a: ",int(gaggvalida),"gmult: ",int(gmult),"xmult: ",gaxmult,temp1a,temp1b,temp1c,temp1d,temp2a,temp2b,temp2c,temp2d)
                del xid
                del xevalid


        return [ID, Ge]



def nsclpxi16evread(fpr):
	'''
	NSCL/FRIB Pixie16 Data Format
	C. Wibisono
	02/16 '24
	Usage:
	To parse the .evt file coming from NSCL/FRIB DAQ. This function parser would return pxi16data object per event.
	Function Argument(s):
	fpr: file pointer object
	Return value: either -1, 0 or pxi16 data type
	-1: end of file
	0: skip for non-physics ring
	pxi16dataobj[sevtmult]: list of pxi16data object. The dimension would denote number of hits for each event.
	'''
	p=Struct("@i")
	q=Struct("@I")
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
			break
		temporary=temporary-(48+2*bsize)
		buff1=fpr.read(4)
		if buff1 == b'':
			break
		buff1int,=q.unpack(buff1)
		chn =buff1int & 0xF
		sln=(buff1int & 0xF0) >> 4
		crn=(buff1int & 0XF00) >> 8
		hlen=(buff1int & 0x1F000) >> 12
		elen=(buff1int & 0x7FFE0000) >> 17
		fcode=(buff1int & 0x80000000) >> 31
		buff2=fpr.read(4)
		buff3=fpr.read(4)
		buff2int,=q.unpack(buff2)
		buff3int,=q.unpack(buff3)
		time=((buff3int & 0xFFFF) << 32) + buff2int			
		#ctime=(buff3int & 0x7FFF0000) >> 16
		cfd_forced= (buff3int & 0x80000000) >> 31
		cfd_source= (buff3int & 0x40000000) >> 30
		cfdtimetemp = (buff3int & 0x3FFF0000) >> 16
		cfdtimetemp2 = cfdtimetemp & 0x3FFF
		#cfdtimetemp2 = cfdtimetemp & 0x7FFF
		if cfd_forced == 0:
			if cfd_source == 1:
				#time = time*2 - 1+ cfdtimetemp2/(2*16384)
				ctime = cfdtimetemp2/16384 - 1
			else:
				#time = time*2 + cfdtimetemp2/(2*16384)
				ctime = cfdtimetemp2/16384 
		else:
			ctime = 0
			#time = time
		
		ctimef= (buff3int & 0x80000000) >> 31
		buff4=fpr.read(4)
		buff4int,=q.unpack(buff4)
		tempenergy= buff4int & 0xFFFF
		tempen=float(tempenergy/2.)
		energy=int(tempen)
		trlen=(buff4int & 0x7FFF0000) >> 16
		trwlen=int(trlen/2.)  
		extra=(buff4int & 0x80000000) >> 31
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
				temp,=q.unpack(buff)
				esum.append(temp)
			if hlen == 16:
				for j in range(0,8,1):
					buff=fpr.read(4)
					temp,=q.unpack(buff)	
					qsum.append(temp)

		if hlen == 12:
			for i in range(0,8,1):
				buff=fpr.read(4)
				temp,=q.unpack(buff)
				qsum.append(temp)

			
		tr=[]
		for i in range(hlen,elen,1):
			buff=fpr.read(4)
			buffint,=q.unpack(buff)
				
			if trwlen !=0:
				temp1=buffint & 0x3FFF
				temp2=(buffint >> 16) & 0x3FFF
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
			temp=nsclpxi16evread(f)
			if temp == -1:
				break
			elif temp == 0:
				continue
			else:
				sevtmult=len(temp) #event multiplicity
				for i in range(0,sevtmult,1):
					print("sevtmult:",sevtmult,"crn:",temp[i].crn,"sln:",temp[i].sln,"chn:",temp[i].chn,"id:",temp[i].iddet,"hlen:",temp[i].hlen,"elen:",temp[i].elen,"trlen:",temp[i].trlen)
					if temp[i].energy >=0 and temp[i].energy < 8192:
						idenmat[temp[i].iddet][temp[i].energy]=idenmat[temp[i].iddet][temp[i].energy]+1
					
	#Write matrix into a file
	matwrite(matfilename,dimy=416,dimx=8192,arr=idenmat,overwrite=overwrite)
