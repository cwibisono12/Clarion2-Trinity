//Convert a 2D Matrix into ROOT 2D Histogram format
//C. Wibisono
//06/02/'23

#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include "TH2.h"
#include "TObject.h"


int main(int argc, char* argv[]){
FILE *fin;
fin=fopen(argv[1],"r");
int dimy=atoi(argv[2]);
int dimx=atoi(argv[3]);
char *p;
strcpy(p,argv[5]);

TFile *f=new TFile(argv[4],"RECREATE");
TH2I *root2D= new TH2I(p,p,dimx,0,dimx,dimy,0,dimy);

int *matin[dimy];

int i,j;
for(i=0;i<dimy;i++){
matin[i]=(int*)malloc(sizeof(int)*dimx);
}

for(i=0;i<dimy;i++){
fread(matin[i],sizeof(int)*dimx,1,fin);
}

for(i=0;i<dimx;i++){
	for(j=0;j<dimy;j++){
		root2D->SetBinContent(i,j,matin[j][i]);	
		}
}
printf("hello\n");

f->Write();

for(i=0;i<dimy;i++){
free(matin[i]);
}

//delete root2D;



return 0;
}
