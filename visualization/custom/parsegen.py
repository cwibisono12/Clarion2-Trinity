#!/usr/bin/env python3
from struct import *
import numpy as np

def matfile(filename,*,dimx=4096,dimy=4096):
	'''
	Matrix file format .spn, .spn2, .sec 
	C. Wibisono
	02/14 '24
	Parameter(s):
	Filename: file pointer object
	dimx = x dimension
	dimy = y dimension
	
	Return:
	arr: ndimensional array with shape (dimy,dimx)
	Usage:
	To parse the matrix file.
	'''
	
	p=Struct("@i")
	arr=np.ndarray(shape=(dimy,dimx),dtype=np.int32)
	for i in range(dimy):
		for j in range(dimx):
			arr[i][j] = 0
	with open(filename,mode='rb') as f:
		for i in range(0,dimy,1):
			for j in range(0,dimx,1):
				buff=f.read(4)
				temp,=p.unpack(buff)
				arr[i][j] = temp

	return arr

def matfile8(filename,*,dimx=4096,dimy=4096):
	'''
	Matrix file for long long int (8 bytes)
	C. Wibisono
	05/17 '24
	Parameter(s):
	Filename: file pointer object
	dimx = x dimension
	dimy = y dimension
	
	Return:
	arr: ndimensional array with shape (dimy,dimx)
	Usage:
	To parse the matrix file.
	'''
	
	p=Struct("@q")
	arr=np.ndarray(shape=(dimy,dimx),dtype=np.int64)
	for i in range(dimy):
		for j in range(dimx):
			arr[i][j] = 0
	with open(filename,mode='rb') as f:
		for i in range(0,dimy,1):
			for j in range(0,dimx,1):
				buff=f.read(8)
				temp,=p.unpack(buff)
				arr[i][j] = temp

	return arr

def matwrite(filename,*,dimy,dimx,arr,overwrite):
	'''
	Matrix writer
	C. Wibisono
	02/18 '24
	Usage:
	To write and or update a matrix into a file.
	Function Arguments
	Filename : file to write/update
	dimy: int, y dimension
	dimx: int, x dimension
	arr: int[dimy][dimx], two dimensional array
	overwrite: int : 1 (to overwrite) or 0 (to append)
	'''
	if overwrite == 1:
		with open(filename,mode='wb') as f:
			for i in range(0,dimy,1):
				for j in range(0,dimx,1):
					temp=arr[i][j]
					f.write(pack("@i",temp))

			print("Completed\n")
	else:
		with open(filename,mode='rb') as f:
			for i in range(0,dimy,1):
				for j in range(0,dimx,1):
					buff=f.read(4)
					temp,=unpack("@i",buff)
					arr[i][j]=arr[i][j]+temp
			print("Complete updating the matrix")
			print("Writing updated matrix")
		with open(filename,mode='wb') as f:
			for i in range(0,dimy,1):
				for j in range(0,dimx,1):
					temp=arr[i][j]
					f.write(pack("@i",temp))

			print("Completed\n")
			



def ev5file(filename):
	'''
	.ev5 file format
	C. Wibisono
	02/14 '24
	Usage:
	To parse the .ev5 file.
	'''
	p=Struct("@b")
	q=Struct("@h")
	with open(filename,mode='rb') as f:
		while(1):
		#for i in range(0,10,1):
			buff1=f.read(1)
			if buff1 == b'':
				break
			gaggvalidp,=p.unpack(buff1)
			buff2=f.read(1)
			if buff2 == b'':
				break
			gaggvalida,=p.unpack(buff2)
			buff3=f.read(1)
			if buff3 == b'':
				break
			gmult,=p.unpack(buff3)

			if gaggvalidp != 0 or gaggvalida != 0:
				buff4=f.read(2)
				if buff4 == b'':
					break
				Ex,=q.unpack(buff4)
			else:
				pass

			for i in range(0,gmult,1):
				xid=[]
				xevalid=[]
				gid,=q.unpack(f.read(2))
				gaold,=q.unpack(f.read(2))
				ganew,=q.unpack(f.read(2))
				gae,=q.unpack(f.read(2))
				gaetrue,=q.unpack(f.read(2))
				gavalidp,=p.unpack(f.read(1))
				gavalidnp,=p.unpack(f.read(1))
				gatime,=q.unpack(f.read(2))
				gaxmult,=p.unpack(f.read(1))
				temp1a,=p.unpack(f.read(1))
				temp1b,=p.unpack(f.read(1))
				temp1c,=p.unpack(f.read(1))
				temp1d,=p.unpack(f.read(1))
				f.read(1)
				temp2a,=q.unpack(f.read(2))
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
				
				print("p: ",int(gaggvalidp),"a: ",int(gaggvalida),"gmult: ",int(gmult),"xmult: ",gaxmult,temp1a,temp1b,temp1c,temp1d,temp2a,temp2b,temp2c,temp2d)
				del xid
				del xevalid
					
					


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
			buff1int,=unpack("@I",buff1)
			chn =buff1int & 0xF
			sln=(buff1int & 0xF0) >> 4
			crn=(buff1int & 0XF00) >> 8
			hlen=(buff1int & 0x1F000) >> 12
			elen=(buff1int & 0x7FFE0000) >> 17
			fcode=(buff1int & 0x80000000) >> 31
			buff2=f.read(4)
			buff3=f.read(4)
			buff2int,=unpack("@I",buff2)
			buff3int,=unpack("@I",buff3)
			temp=(buff3int & 0xFFFF) << 32
			time= temp + buff2int		
			ctime=(buff3int & 0x7FFF0000) >> 16
			ctimef=(buff3int & 0x80000000) >> 31
			buff4=f.read(4)
			buff4int,=unpack("@I",buff4)
			energy=buff4int & 0xFFFF
			trlen=(buff4int & 0x7FFF0000) >> 16
			trwlen=int(trlen/2.)  
			extra=(buff4int & 0x80000000) >> 31
			
			if hlen == 4 and trwlen == 0:
				continue
			#f.seek(-16,1)
			esum=[]
			qsum=[]
			if hlen == 8 or hlen == 16:
				for i in range(0,4,1):
					buff=f.read(4)
					esum.append(unpack("@I",buff))
				if hlen == 16:
					for j in range(0,8,1):
						buff=f.read(4)	
						qsum.append(unpack("@I",buff))

			if hlen == 12:
				for i in range(0,8,1):
					buff=f.read(4)
					qsum.append(unpack("@I",buff))

			
			tr=[]
			for i in range(hlen,elen,1):
				buff=f.read(4)
				buffint,=unpack("@I",buff)				
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

			print("chn:",chn,"sln:",sln,"crn:",crn,"hlen:",hlen,"elen:",elen,"time:",time,"buff2:",buff2int,"temp:",temp)


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
				buff1int,=unpack("@I",buff1)
				chn =buff1int & 0xF
				sln=(buff1int & 0xF0) >> 4
				crn=(buff1int & 0XF00) >> 8
				hlen=(buff1int & 0x1F000) >> 12
				elen=(buff1int & 0x7FFE0000) >> 17
				fcode=(buff1int & 0x80000000) >> 31
				buff2=f.read(4)
				buff3=f.read(4)
				buff2int,=unpack("@I",buff2)
				buff3int,=unpack("@I",buff3)
				time=((buff3int & 0xFFFF) << 32) + buff2int			
				ctime=(buff3int & 0x7FFF0000) >> 16
				ctimef=(buff3int & 0x80000000) >> 31
				buff4=f.read(4)
				buff4int,=unpack("@I",buff4)
				energy=buff4int & 0xFFFF
				trlen=(buff4int & 0x7FFF0000) >> 16
				trwlen=int(trlen/2.)  
				extra=(buff4int & 0x80000000) >> 31
			
				if hlen == 4 and trwlen == 0:
					continue
				#f.seek(-16,1)
				esum=[]
				qsum=[]
				if hlen == 8 or hlen == 16:
					for i in range(0,4,1):
						buff=f.read(4)
						esum.append(unpack("@I",buff))
				if hlen == 16:
					for j in range(0,8,1):
						buff=f.read(4)	
						qsum.append(unpack("@I",buff))

				if hlen == 12:
					for i in range(0,8,1):
						buff=f.read(4)
						qsum.append(unpack("@I",buff))

			
				tr=[]
				for i in range(hlen,elen,1):
					buff=f.read(4)
					buffint,=unpack("@I",buff)				
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
				print("chn:",chn,"sln:",sln,"crn:",crn,"hlen:",hlen,"elen:",elen,"energy:",energy)



if __name__ == "__main__":
	'''
	Below is an example how to use this module.
	The example below demonstrates how to read the raw file coming from NSCL/FRIB DAQ.
	'''
	import sys
	filename=sys.argv[1]
	nsclpxi16file(filename)
