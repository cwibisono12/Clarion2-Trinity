#include <stdio.h>
#include <stdlib.h>

//Convert 2D Matrices into another 2D matrices with different dimensions:
//C. Wibisono
//04/19/2023

int main(int argc, char* argv[]){

FILE *fin, *fout;
fin=fopen(argv[1],"r");
fout=fopen(argv[2],"w");
int sizeY=atoi(argv[3]);
int sizeX=atoi(argv[4]);

int Xfout=atoi(argv[5]);
int Yfout=atoi(argv[6]);

int *matin[sizeY];
int *matout[sizeX];

int i,j;

//Allocate the Memory:
for(i=0;i<sizeY;i++){
matin[i]=(int*)malloc(sizeof(int)*sizeX);
}

//Reading Matrices Input File
for(j=0;j<sizeY;j++){
fread(matin[j],sizeof(int)*sizeX,1,fin);
}

for(j=0;j<Yfout;j++){
matout[j]=(int*)malloc(sizeof(int)*Xfout);
}

for(i=0;i<Yfout;i++){
	for(j=0;j<Xfout;j++){
	matout[i][j]=matin[i][j];
		}
}

//Deallocate Matrices Input File:
for(j=0;j<sizeY;j++){
	free(matin[j]);
}

//Write Output Matrices into a file:
for(i=0;i<Yfout;i++){
fwrite(matout[i],sizeof(int)*Xfout,1,fout);
}

//Deallocate Matrices Output File:
for(j=0;j<Yfout;j++){
free(matout[j]);
}
return 0;
}
