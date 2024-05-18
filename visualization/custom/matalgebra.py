#!/usr/bin/env python3

from parsegen import matfile as mat
from parsegen import matfile8 as mat2
import numpy as np

def addmat4(filemat1, filemat2, * ,dimy=4096,dimx=4096, scalefactor = 1):	
	'''
	Function to Add 2 Matrix Files for integer 4 bytes 
	C. Wibisono
	05/17 '24
	Parameter(s):
	filemat1  = matrix file1 pointer
	filemat2 = matrix file2 pointer
	dimy = y-dimension
	dimx = x-dimension
	scalefactor = scaling factor
	Return:
	arr[dimy][dimx]
	'''
	m1 = mat(filemat1, dimx=dimx, dimy=dimy)
	m2 = mat(filemat2, dimx=dimx, dimy=dimy)
	m3 = np.ndarray(shape=(dimy,dimx),dtype = np.int32)

	for i in range(dimy):
		for j in range(dimx):
			m3[i][j]=0

	for k1 in range(dimy):
		for k2 in range(dimx):
			m3[k1][k2] = m1[k1][k2] + m2[k1][k2]*scalefactor

	return m3

def submat4(filemat1, filemat2, * ,dimy=4096,dimx=4096, scalefactor = 1):
	'''
	Function to Subtract 2 Matrix Files for integer 4 bytes 
	C. Wibisono
	05/17 '24
	Parameter(s):
	filemat1  = matrix file1 pointer
	filemat2 = matrix file2 pointer
	dimy = y-dimension
	dimx = x-dimension
	scalefactor = scaling factor
	Return:
	arr[dimy][dimx]
	'''
	m1 = mat(filemat1, dimx=dimx, dimy=dimy)
	m2 = mat(filemat2, dimx=dimx, dimy=dimy)
	m3 = np.ndarray(shape=(dimy,dimx),dtype = np.int32)
	
	for i in range(dimy):
		for j in range(dimx):
			m3[i][j]=0
	
	for k1 in range(dimy):
		for k2 in range(dimx):
			m3[k1][k2] = m1[k1][k2] - m2[k1][k2]*scalefactor
	
	return m3

	
def addmat8(filemat1, filemat2, * ,dimy=4096,dimx=4096, scalefactor = 1):
	'''
	Function to Add 2 Matrix Files for integer 8 bytes 
	C. Wibisono
	05/17 '24
	Parameter(s):
	filemat1  = matrix file1 pointer
	filemat2 = matrix file2 pointer
	dimy = y-dimension
	dimx = x-dimension
	scalefactor = scaling factor
	Return:
	arr[dimy][dimx]
	'''
	m1 = mat2(filemat1, dimx=dimx, dimy=dimy)
	m2 = mat2(filemat2, dimx=dimx, dimy=dimy)
	m3 = np.ndarray(shape=(dimy,dimx),dtype = np.int64)
	
	for i in range(dimy):
		for j in range(dimx):
			m3[i][j]=0
	
	for k1 in range(dimy):
		for k2 in range(dimx):
			m3[k1][k2] = m1[k1][k2] + m2[k1][k2]*scalefactor
	
	return m3

def submat8(filemat1, filemat2, * ,dimy=4096,dimx=4096, scalefactor = 1):
	'''
	Function to Subtract 2 Matrix Files for integer 8 bytes 
	C. Wibisono
	05/17 '24
	Parameter(s):
	filemat1  = matrix file1 pointer
	filemat2 = matrix file2 pointer
	dimy = y-dimension
	dimx = x-dimension
	scalefactor = scaling factor
	Return:
	arr[dimy][dimx]
	'''
	m1 = mat2(filemat1, dimx=dimx, dimy=dimy)
	m2 = mat2(filemat2, dimx=dimx, dimy=dimy)
	m3 = np.ndarray(shape=(dimy,dimx),dtype = np.int64)
	
	for i in range(dimy):
		for j in range(dimx):
			m3[i][j]=0
	
	for k1 in range(dimy):
		for k2 in range(dimx):
			m3[k1][k2] = m1[k1][k2] - m2[k1][k2]*scalefactor
	
	return m3

if __name__ =="__main__":
	import sys
	from parsegen import matwrite 

	f1 = sys.argv[1] #matfile1
	f2 = sys.argv[2] #matfile2
	fout = sys.argv[3] #matfileoutput
	scale = int(sys.argv[4]) #scaling factor
	dimy = int(sys.argv[5]) #y-dimension
	dimx = int(sys.argv[6]) #x-dimension

	#add f1 + f2:
	matrix = addmat4(f1,f2,dimy=dimy,dimx=dimx,scalefactor=scale)
	
	#write matrix into fout:
	matwrite(fout,dimy=dimy,dimx=dimx,arr=matrix,overwrite = 1)
