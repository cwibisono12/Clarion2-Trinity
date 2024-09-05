#include <stdio.h>
#include <stdlib.h>
#include "global.h"

//Driver Code for pxi16filter.c

int main(int argc, char* argv[]){

FILE *fin, *fout;
fin = fopen(argv[1],"r"); //Original Pxi16 Time-Ordered List Mode File
fout = fopen(argv[2],"w"); //Filtered Pxi16 List Mode Data

struct pxi16event pxi16evt[MAX_ID]={0};
struct subevent subevt[MAX_ID]={0}; 
unsigned int sub[MAX_SUB_LENGTH];

int sevtmult;
int i,j,k;

while(1){
int Ge = 0;
int GAGG = 0;

sevtmult = tofilter(sub,subevt,fin,pxi16evt); 
	for(i=0;i<sevtmult;i++){
		if (pxi16evt[i].iddet >=0 && pxi16evt[i].iddet <=57) Ge = Ge+1; 
		if (pxi16evt[i].iddet >=64 && pxi16evt[i].iddet <=191) GAGG = GAGG + 1;
		}
	
	//Only Allow event that contains gamma hits to be into the filtered pxi16 list mode data:
	if(Ge > 0){
	for(i=0;i<sevtmult;i++){
		j = pxi16evt[i].elen;
		for(k=0;k<j;k++){
		fwrite(&pxi16evt[i].arr[k], sizeof(int), 1, fout);
		}
	}


	}	
if(sevtmult == 0) break;
}

fclose(fin);
fclose(fout);

return 0;
}
