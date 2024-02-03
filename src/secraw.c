#include "global.h"
#include "evtnscl.h"

void secraw(int len, struct subevent *subevt){

int i, flagim=0, crosscint=0, crossmsx40=0;
int indim,indcrosscint,indcrossmsx40;

for(i=0;i <len;i++){
if((subevt[i].id >=0 && subevt[i].id <=416) &&(subevt[i].energy >=0 && subevt[i].energy<=8192)){
rawhisto[subevt[i].id][subevt[i].energy]++;
}

if(subevt[i].crn ==0){
hitcr1[subevt[i].sln][subevt[i].chn]++;
}

if(subevt[i].crn==1){
hitcr2[subevt[i].sln][subevt[i].chn]++;
}


//Time of Flight:
//Image Scint:
if(subevt[i].id == 26 || subevt[i].id == 27) {flagim++;
indim=i;
}
//Cross Scint:
if(subevt[i].id == 40 || subevt[i].id == 24) {crosscint++;
indcrosscint=i;
}

//For energy loss through upstream PIN Si detector before the Implant Detector:
//if(subevt[i].id == 25) {crossmsx40++;
if(subevt[i].id == 32 || subevt[i].id == 33) {crossmsx40++;
indcrossmsx40=i;
//printf("energy: %d\n",subevt[i].energy);
	}


}

int timediff=0;
int timediffnew=0;
//if((flagim==1) && (crosscint == 1 && crossmsx40 >0))
//printf("flagim: %d crossscint: %d crossmsx40: %d dE: %d\n",flagim,crosscint,crossmsx40,subevt[indcrossmsx40].energy);
//long long int tdiff1;
int temp;

if((flagim == 1 && crosscint == 1) &&  crossmsx40 > 0){
//temp=(subevt[indiml].time+subevt[indimr].time)/2.;
//printf("flagiml: %d,flagimr: %d, crossscint: %d crossmsx40: %d\n",flagiml,flagimr,crosscint,crossmsx40);

//tdiff1=(long long int) temp;


if(flagim==1){
timediff=(subevt[indcrosscint].time-subevt[indim].time);
temp=timediff*1.;
timediffnew=(int) temp;
}
int energy=subevt[indcrossmsx40].energy;
printf("mult: %d timediffnew: %d index: %d energymsx40: %d\n",len,timediffnew, indcrossmsx40, energy);
if((timediffnew>0 && timediffnew <4096) && (energy > 0 && energy <4096))
implantpid[energy][timediffnew]++;
}


}
