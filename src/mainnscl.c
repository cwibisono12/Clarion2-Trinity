#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "evtnscl.h"
#include "evtnscl_dec.h"


int main(int argc, char* argv[]){
FILE *fp;
fp=fopen(argv[1],"r");
int overwrite=0;
int option=atoi(argv[2]);
int ab,bc;
for(ab=0;ab<416;ab++){
if((rawhisto[ab]=(int*)malloc(8192*sizeof(int)))==NULL){
printf("\nError,memory not allocated");
exit(1);
}
}

for(bc=0;bc<20;bc++){
if((hitcr1[bc]=(int*)malloc(16*sizeof(int)))==NULL){
printf("\nError,memory not allocated");
exit(1);
}
if((hitcr2[bc]=(int*)malloc(16*sizeof(int)))==NULL){
printf("\nError,memory not allocated");
exit(1);
}
}
/*
int ab;
for(ab=0;ab<4096;ab++){
	if((pid_qdc21[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc22[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc23[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc24[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc25[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc26[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc27[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc28[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc29[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc210[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc41[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc42[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc43[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc44[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc45[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc46[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc47[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc48[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc49[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc410[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc411[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc412[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc413[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc414[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc415[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
	if((pid_qdc416[ab]=(int*)malloc(4096*sizeof(int)))== NULL){
	printf("\nError,memory not allocated");
	exit(1);
	}
}

*/


unsigned int sub[MAX_SUB_LENGTH];
memset(sub,0,sizeof(sub));
struct subevent subevt[1000];
evtreader(sub,subevt,option,fp);
write_data4dyn(RAW_HISTO,rawhisto,8192,416,overwrite);
write_data4dyn(HITCR1,hitcr1,16,20,overwrite);
write_data4dyn(HITCR2,hitcr2,16,20,overwrite);
/*
write_data4dyn(PID_QDC21, pid_qdc21, 4096, 4096, overwrite);
write_data4dyn(PID_QDC22, pid_qdc22, 4096, 4096, overwrite);
write_data4dyn(PID_QDC23, pid_qdc23, 4096, 4096, overwrite);
write_data4dyn(PID_QDC24, pid_qdc24, 4096, 4096, overwrite);
write_data4dyn(PID_QDC25, pid_qdc25, 4096, 4096, overwrite);
write_data4dyn(PID_QDC26, pid_qdc26, 4096, 4096, overwrite);
write_data4dyn(PID_QDC27, pid_qdc27, 4096, 4096, overwrite);
write_data4dyn(PID_QDC28, pid_qdc28, 4096, 4096, overwrite);
write_data4dyn(PID_QDC29, pid_qdc29, 4096, 4096, overwrite);
write_data4dyn(PID_QDC210, pid_qdc210, 4096, 4096, overwrite);
write_data4dyn(PID_QDC41, pid_qdc41, 4096, 4096, overwrite);
write_data4dyn(PID_QDC42, pid_qdc42, 4096, 4096, overwrite);
write_data4dyn(PID_QDC43, pid_qdc43, 4096, 4096, overwrite);
write_data4dyn(PID_QDC44, pid_qdc44, 4096, 4096, overwrite);
write_data4dyn(PID_QDC45, pid_qdc45, 4096, 4096, overwrite);
write_data4dyn(PID_QDC46, pid_qdc46, 4096, 4096, overwrite);
write_data4dyn(PID_QDC47, pid_qdc47, 4096, 4096, overwrite);
write_data4dyn(PID_QDC48, pid_qdc48, 4096, 4096, overwrite);
write_data4dyn(PID_QDC49, pid_qdc49, 4096, 4096, overwrite);
write_data4dyn(PID_QDC410, pid_qdc410, 4096, 4096, overwrite);
write_data4dyn(PID_QDC411, pid_qdc411, 4096, 4096, overwrite);
write_data4dyn(PID_QDC412, pid_qdc412, 4096, 4096, overwrite);
write_data4dyn(PID_QDC413, pid_qdc413, 4096, 4096, overwrite);
write_data4dyn(PID_QDC414, pid_qdc414, 4096, 4096, overwrite);
write_data4dyn(PID_QDC415, pid_qdc415, 4096, 4096, overwrite);
write_data4dyn(PID_QDC416, pid_qdc416, 4096, 4096, overwrite);
*/
return 0;
}
