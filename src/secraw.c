#include "global.h"
#include "evtnscl.h"

void secraw(int len, struct subevent *subevt){

int i;
for(i=0;i <len;i++){
if((subevt[i].id >=0 && subevt[i].id <=416) &&(subevt[i].energy >=0 && subevt[i].energy<=8192)){
rawhisto[subevt[i].id][subevt[i].energy]++;
}

if(subevt[i].crn ==0){
hitcr1[subevt->sln][subevt->chn]++;
}

if(subevt[i].crn==1){
hitcr2[subevt->sln][subevt->chn]++;
}

}
}
