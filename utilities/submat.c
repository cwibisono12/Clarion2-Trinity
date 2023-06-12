#include <stdio.h>
#include <stdlib.h>

//Subtract two Matrices with Equal Dimension
//C. Wibisono
//06/12 '23
//How to Use:
//./submat gg_prompt.spn2 gg_nprompt.spn2 gg_prompt_clean.spn2 5000 5000


int main(int argc, char* argv[]){
FILE *fin1,*fin2, *fout;
fin1=fopen(argv[1],"r");
fin2=fopen(argv[2],"r");
fout=fopen(argv[3],"w");

int dimy=atoi(argv[4]);
int dimx=atoi(argv[5]);

int *matin1[dimy];
int *matin2[dimy];
int *matout[dimy];

int i,j;
for(i=0;i<dimy;i++){
matin1[i]=(int*)malloc(sizeof(int)*dimx);
matin2[i]=(int*)malloc(sizeof(int)*dimx);
matout[i]=(int*)malloc(sizeof(int)*dimx);
}

for(i=0;i<dimy;i++){
fread(matin1[i],sizeof(int)*dimx,1,fin1);
fread(matin2[i],sizeof(int)*dimx,1,fin2);
}

for(i=0;i<dimy;i++){
	for(j=0;j<dimx;j++){
		matout[i][j]=matin1[i][j]-matin2[i][j];
		}	
}

for(i=0;i<dimy;i++){
fwrite(matout[i],sizeof(int)*dimx,1,fout);
}

for(i=0;i<dimy;i++){
free(matin1[i]);
free(matin2[i]);
free(matout[i]);
}

fclose(fin1);
fclose(fin2);
fclose(fout);
return 0;
}
