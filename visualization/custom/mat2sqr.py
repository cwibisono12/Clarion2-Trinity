#!/usr/bin/env python3

from parsegen import matfile, matwrite

from struct import *
import numpy as np


def writesqr(fin,fout,*,dimy=5000,dimx=5000):
	'''
	Usage:
	To convert the original matrix file format into local square matrix file format (Internal Use only).
	Parameter(s):
	fin: original file format
	fout: local file format
	dimy: y dimension
	dimx: x dimension
	'''

	p=Struct("@i")
	q=Struct("@f")
	a=12
	tempbuff=[dimx,0,0,0]
	dbuff=4*dimy*dimx

	#read original matrix file:
	print("Loading the input file\n")
	arr=matfile(fin,dimx=dimx,dimy=dimy)
	print("Complete Loading the input file\n")
	#write local square matrix file:
	with open(fout, mode='wb') as f:
		f.write(p.pack(a))
		f.write(p.pack(tempbuff[0]))
		f.write(p.pack(tempbuff[1]))
		f.write(p.pack(tempbuff[2]))
		f.write(p.pack(a))
		f.write(p.pack(dbuff))
		print("Writing matrix into an output file\n")
		for i in range(dimy):
			for j in range(dimx):
				temp=float(arr[i][j])
				f.write(q.pack(temp))
				#if arr[i][j] !=0:
				#	print("i:",i,"j:",j,"value:",temp)

		f.write(p.pack(dbuff))

	print("Process Completed.")


def readsqr(fin,*,dimy=5000,dimx=5000):
	'''
	Usage:
	To read the local square matrix file
	Parameter(s):
	fin: file input pointer object
	dimy: y-dimension
	dimx: x-dimension
	arr: 2D[dimy][dimx] matrix
	'''
	p=Struct("@i")
	q=Struct("@f")
	a=12
	tempbuff=[dimx,0,0,0]
	dbuff=4*dimy*dimx
	arr=np.ndarray(shape=(dimy,dimx),dtype=np.int32)

	for i in range(dimy):
		for j in range(dimx):
			arr[i][j]=0

	with open(fin, mode='rb') as f:
		f.read(4)
		f.read(4*3)
		f.read(4)
		f.read(4)
		for i in range(dimy):
			for j in range(dimx):
				temp,=q.unpack(f.read(4))
				arr[i][j]=int(temp)
				if arr[i][j] !=0:
					print("i:",i,"j:",j,"value:",temp)

		f.read(4)

	print("Completed.")
	return arr

if __name__ == "__main__":
	import sys
	filein = sys.argv[1]
	fileout = sys.argv[2]
	writesqr(filein,fileout,dimy=5000,dimx=5000)
