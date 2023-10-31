//Unpack Pixie16 Digitizer

#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "spec_dec.h"
#include "histoevt.h"
#include "histoevt_dec.h"




void evtreader(unsigned int sub[], struct subevent *subevt, FILE *fpr){
      /////////////////////////////////
      // UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////
      int i=0,j=0,k=0;
        
       while(1){ 
        //read 4-byte header
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        subevt->chn = sub[0] & 0xF;
        subevt->sln = (sub[0] & 0xF0) >> 4;
        subevt->crn = (sub[0] & 0xF00) >> 8;
        subevt->id = subevt->crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (subevt->sln - BOARD_START)*MAX_CHANNELS_PER_BOARD + subevt->chn;   
        subevt->hlen = (sub[0] & 0x1F000) >> 12;
        subevt->elen = (sub[0] & 0x7FFE0000) >> 17;
        subevt->fcode = (sub[0] & 0x80000000) >> 31;
        subevt->time = ( (long long int)(sub[2] & 0xFFFF) << 32) + sub[1];
        subevt->ctime = (sub[2] & 0x7FFF0000) >> 16;
        subevt->ctimef = (sub[2] & 0x80000000) >> 31;
        subevt->energy = (sub[3] & 0xFFFF);
        subevt->trlen = (sub[3] & 0x7FFF0000) >> 16;
        subevt->trwlen = subevt->trlen / 2;
        subevt->extra = (sub[3] & 0x80000000) >> 31;       
 
         
        //continue on if no trace, esum, or qsum
        if (subevt->hlen==HEADER_LENGTH && subevt->trwlen==0 ) {
            continue;
        }
        //more data than just the header; read entire sub event
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        if (fread(sub, sizeof(int)*subevt->elen, 1, fpr) != 1) break;
                              
        //trace
        k=0;
        for (i = subevt->hlen; i < subevt->elen; i++) {      
            subevt->tr[i - subevt->hlen + k] = sub[i] & 0x3FFF; // the upper 2 bits/16 bits are filled with 0s
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0xFFFF; //C. W  
           subevt->tr[i - subevt->hlen + k + 1] = (sub[i]>>16) & 0x3FFF;
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0xFFFF; //C. W
            k=k+1;
        } 


       	//printf("id: %d hlen: %d elen: %d\n",subevt->id,subevt->hlen,subevt->elen); 
	if(subevt->trlen != 0){
	int l,m;
	for(m=0;m<8;m++){
	subevt->qsum[m]=0;
	}
	for(l=0;l<31;l++){
	subevt->qsum[0]=subevt->qsum[0]+subevt->tr[l];
	}
	for(l=31;l<60;l++){
	subevt->qsum[1]=subevt->qsum[1]+subevt->tr[l];
	}
	for(l=60;l<75;l++){
	subevt->qsum[2]=subevt->qsum[2]+subevt->tr[l];
	}
	for(l=75;l<95;l++){
	subevt->qsum[3]=subevt->qsum[3]+subevt->tr[l];
	}
	for(l=95;l<105;l++){
	subevt->qsum[4]=subevt->qsum[4]+subevt->tr[l];
	}
	for(l=105;l<160;l++){
	subevt->qsum[5]=subevt->qsum[5]+subevt->tr[l];
	}
	for(l=160;l<175;l++){
	subevt->qsum[6]=subevt->qsum[6]+subevt->tr[l];
	}
	for(l=175;l<200;l++){
	subevt->qsum[7]=subevt->qsum[7]+subevt->tr[l];
	}
	}

        
     // if (subevt[sevtmult].id == 4 && subevt[sevtmult].fcode == 1) DB(subevt[sevtmult].tr);
            
        //continue if no esum or qsum   
        if (subevt->hlen==HEADER_LENGTH) {
	    pidevt(subevt);
            continue;
        }
        
        //esum
        if (subevt->hlen==8 || subevt->hlen==16) { 
            for (i=4; i < 8; i++) {
                subevt->esum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt->hlen==12) { 
            for (i=4; i < 12; i++) {
                subevt->qsum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt->hlen==16) { 
            for (i=8; i < 16; i++) {
                subevt->qsum[i-8] = sub[i];
            }
        }    
     
	//if(subevt->id >=60 && subevt-> id <=120)
	//printf("id: %d hlen: %d elen: %d q1: %d q2: %d q3: %d\n",subevt->id,subevt->hlen,subevt->elen,subevt->qsum[0],subevt->qsum[1],subevt->qsum[2]);
	pidevt(subevt);   
    }
 

}



