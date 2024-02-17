#!/usr/bin/env python3
from struct import *

class gdetectorsum:
	def __init__(self,detid,angleold,anglenew,energy,etrue,validp,validnp,time,xmult,xid,xevalid):
		self.detid=detid
		self.angleold=angleold
		self.anglenew=anglenew
		self.energy=energy
		self.etrue=etrue
		self.validp=validp
		self.validnp=validnp
		self.time=time
		self.xmult=xmult
		self.xid=xid
		self.xevalid=xevalid
	

def matfile(filename,*,dimx=4096,dimy=4096):
	'''
	Matrix file format .spn, .spn2, .sec 
	C. Wibisono
	02/14 '24
	Usage:
	To parse the matrix file.
	'''
	with open(filename,mode='rb') as f:
		for i in range(0,dimy,1):
			for j in range(0,dimx,1):
				buff=f.read(4)
				temp,=unpack("@i",buff)
				if temp != 0:
					print("i: ",i,"j: ",j,"val: ",temp)



def ev5file(filename):
	'''
	.ev5 file format
	C. Wibisono
	02/14 '24
	Usage:
	To parse the .ev5 file.
	'''
	with open(filename,mode='rb') as f:
		while(1):
		#for i in range(0,10,1):
			buff1=f.read(1)
			if buff1 == b'':
				break
			gaggvalidp,=unpack("@b",buff1)
			buff2=f.read(1)
			if buff2 == b'':
				break
			gaggvalida,=unpack("@b",buff2)
			buff3=f.read(1)
			if buff3 == b'':
				break
			gmult,=unpack("@b",buff3)
			buff4=f.read(2)
			if buff4 == b'':
				break
			Ex,=unpack("@h",buff4)
			ge=[]
			#print("p: ",gaggvalidp,"a: ",gaggvalida,"gmult: ",gmult)
			for i in range(0,gmult,1):
				xid=[]
				xevalid=[]
				gid,=unpack("@h",f.read(2))
				gaold,=unpack("@h",f.read(2))
				ganew,=unpack("@h",f.read(2))
				gae,=unpack("@h",f.read(2))
				gaetrue,=unpack("@h",f.read(2))
				gavalidp,=unpack("@b",f.read(1))
				gavalidnp,=unpack("@b",f.read(1))
				gatime,=unpack("@h",f.read(2))
				gaxmult,=unpack("@b",f.read(1))
				temp1a,=unpack("@b",f.read(1))
				temp1b,=unpack("@b",f.read(1))
				temp1c,=unpack("@b",f.read(1))
				temp1d,=unpack("@b",f.read(1))
				f.read(1)
				temp2a,=unpack("@h",f.read(2))
				temp2b,=unpack("@h",f.read(2))
				temp2c,=unpack("@h",f.read(2))
				temp2d,=unpack("@h",f.read(2))
				xid.append(temp1a)
				xid.append(temp1b)
				xid.append(temp1c)
				xid.append(temp1d)
				xevalid.append(temp2a)
				xevalid.append(temp2b)
				xevalid.append(temp2c)
				xevalid.append(temp2d)
				
				ge.append(gdetectorsum(gid,gaold,ganew,gae,gaetrue,gavalidp,gavalidnp,gatime,gaxmult,xid,xevalid))
				print("p: ",int(gaggvalidp),"a: ",int(gaggvalida),"gmult: ",int(gmult),"xmult: ",gaxmult,temp1a,temp1b,temp1c,temp1d,temp2a,temp2b,temp2c,temp2d)
				del xid
				del xevalid
					
					
			del ge


def pxi16file(filename):
	'''
	Original Pixie16 Data Format
	C. Wibisono
	02/16 '24	
	Usage:
	To parse the .evt file coming from original pixie16 data structure
	'''
	with open(filename,mode='rb') as f:
		while(1):	
			buff1=f.read(4)
			if buff1 == b'':
				break
			buff1int,=unpack("@i",buff1)
			chn =int(buff1int & 0xF)
			sln=int((buff1int & 0xF0) >> 4)
			crn=int((buff1int & 0XF00) >> 8)
			hlen=int((buff1int & 0x1F000) >> 12)
			elen=int((buff1int & 0x7FFE0000) >> 17)
			fcode=int((buff1int & 0x80000000) >> 31)
			buff2=f.read(4)
			buff3=f.read(4)
			buff2int,=unpack("@i",buff2)
			buff3int,=unpack("@i",buff3)
			time=int(((buff3int & 0xFFFF) << 32) + buff2int)			
			ctime=int((buff3int & 0x7FFF0000) >> 16)
			ctimef=int((buff3int & 0x80000000) >> 31)
			buff4=f.read(4)
			buff4int,=unpack("@i",buff4)
			energy=int (buff4int & 0xFFFF)
			trlen=int ((buff4int & 0x7FFF0000) >> 16)
			trwlen=int(trlen/2.)  
			extra=int((buff4int & 0x80000000) >> 31)
			
			if hlen == 4 and trwlen == 0:
				continue
			#f.seek(-16,1)
			esum=[]
			qsum=[]
			if hlen == 8 or hlen == 16:
				for i in range(0,4,1):
					buff=f.read(4)
					esum.append(unpack("@i",buff))
				if hlen == 16:
					for j in range(0,8,1):
						buff=f.read(4)	
						qsum.append(unpack("@i",buff))

			if hlen == 12:
				for i in range(0,8,1):
					buff=f.read(4)
					qsum.append(unpack("@i",buff))

			
			tr=[]
			for i in range(hlen,elen,1):
				buff=f.read(4)
				buffint,=unpack("@i",buff)				
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


def nsclpxi16file(filename):
	'''
	NSCL/FRIB Pixie16 Data Format
	C. Wibisono
	02/16 '24
	Usage:
	To parse the .evt file coming from NSCL/FRIB DAQ
	'''
	with open(filename,mode='rb') as f:
		while(1):

			buffring=f.read(4)
			if buffring == b'':
				break
			rihsize,=unpack("@i",buffring)
			bufftype=f.read(4)
			rihtype,=unpack("@i",bufftype)
			if rihtype != 30:
				f.seek(rihsize-8,1)
				continue
			buffsize=f.read(4)
			ribhsize,=unpack("@i",buffsize)
			buffht=f.read(8)
			ribht,=unpack("@q",buffht)
			buffsid=f.read(4)
			ribhsid,=unpack("@i",buffsid)
			buffhbt=f.read(4)
			ribhbt,=unpack("@i",buffhbt)
			buffbhsize=f.read(4)
			ribhbsize,=unpack("@i",buffbhsize)
			temporary=ribhbsize-4
			while(temporary>0):
				f.read(8) #frtmp
				f.read(4) #sid
				f.read(4) #payload
				f.read(4) #bt
				f.read(4) #rihsize
				f.read(4) #rihtype
				f.read(4) #ribhsize
				f.read(8) #temp
				f.read(4) #sid2
				f.read(4) #bt2
				buffbsize=f.read(4) #bsize
				bsize,=unpack("@i",buffbsize)
				buffdevice=f.read(4)
				if buffdevice == b'':
					break
				temporary=temporary-(48+2*bsize)
				buff1=f.read(4)
				if buff1 == b'':
					break
				buff1int,=unpack("@i",buff1)
				chn =int(buff1int & 0xF)
				sln=int((buff1int & 0xF0) >> 4)
				crn=int((buff1int & 0XF00) >> 8)
				hlen=int((buff1int & 0x1F000) >> 12)
				elen=int((buff1int & 0x7FFE0000) >> 17)
				fcode=int((buff1int & 0x80000000) >> 31)
				buff2=f.read(4)
				buff3=f.read(4)
				buff2int,=unpack("@i",buff2)
				buff3int,=unpack("@i",buff3)
				time=int(((buff3int & 0xFFFF) << 32) + buff2int)			
				ctime=int((buff3int & 0x7FFF0000) >> 16)
				ctimef=int((buff3int & 0x80000000) >> 31)
				buff4=f.read(4)
				buff4int,=unpack("@i",buff4)
				energy=int (buff4int & 0xFFFF)
				trlen=int ((buff4int & 0x7FFF0000) >> 16)
				trwlen=int(trlen/2.)  
				extra=int((buff4int & 0x80000000) >> 31)
			
				if hlen == 4 and trwlen == 0:
					continue
				#f.seek(-16,1)
				esum=[]
				qsum=[]
				if hlen == 8 or hlen == 16:
					for i in range(0,4,1):
						buff=f.read(4)
						esum.append(unpack("@i",buff))
				if hlen == 16:
					for j in range(0,8,1):
						buff=f.read(4)	
						qsum.append(unpack("@i",buff))

				if hlen == 12:
					for i in range(0,8,1):
						buff=f.read(4)
						qsum.append(unpack("@i",buff))

			
				tr=[]
				for i in range(hlen,elen,1):
					buff=f.read(4)
					buffint,=unpack("@i",buff)				
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
				print("chn:",chn,"sln:",sln,"crn:",crn,"hlen:",hlen,"elen:",elen,"energy:",energy)



if __name__ == "__main__":
	import sys
	filename=sys.argv[1]
	pxi16file(filename)
	#nsclpxi16file(filename)
	#ev5file(filename)
