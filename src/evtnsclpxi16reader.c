//Unpack Pixie16 Digitizer
//Added Support to include NSCL evt version
//This version of NSCL evt already included event-building
//C. Wibisono 01/27 '24

#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "evtnscl.h"
#include "nsclheader.h"



void evtreader(unsigned int sub[], struct subevent *subevt,int option, FILE *fpr){
      int i=0,j=0,k=0;
      int rihsize=0;
      int rihtype=0;
      int ribhsize=0;
      long long int ribht=0;
      int ribhsid=0;
      int ribhbt=0;
      int ribhbsize=0;
      struct fragment{
	long long int frtmp;
	int sid;
	int payload;
	int bt;
	int rihsize;
	int rihtype;
	int ribhsize;
	long long int tmp;
	int sid2;
	int bt2;
	int bsize;
	int device;
};
	typedef struct fragment FRAG;
	FRAG frag[1000]={0};
	while(1){ 	

	if ( fread(&rihsize,sizeof(int),1,fpr)!=1) break;
	if ( fread(&rihtype,sizeof(int),1,fpr)!=1) break;
	if(rihtype !=30) {fseek(fpr,rihsize-8,SEEK_CUR);
        continue;}
	if( fread(&ribhsize,sizeof(int),1,fpr) !=1) break;
	if( fread(&ribht,sizeof(long long int),1,fpr) !=1) break;
	if( fread(&ribhsid,sizeof(int),1,fpr)!=1) break;
	if( fread(&ribht,sizeof(int),1,fpr)!=1)break;
	if(fread(&ribhbsize,sizeof(int),1,fpr)!=1) break;
	int temporary=ribhbsize-4;
	int ind=0;
	int num=100;
	//iterate for each fragment for a given ring //C.W
	while(temporary>0){
	if(fread(&frag[ind],sizeof(FRAG),1,fpr)!=1) break;
	temporary=temporary-(48+2*frag[ind].bsize);
	
	
	
	//read 4-byte header
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        subevt[ind].chn = sub[0] & 0xF;
        subevt[ind].sln = (sub[0] & 0xF0) >> 4;
        subevt[ind].crn = (sub[0] & 0xF00) >> 8;
        subevt[ind].id = subevt[ind].crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (subevt[ind].sln - BOARD_START)*MAX_CHANNELS_PER_BOARD + subevt[ind].chn;   
        subevt[ind].hlen = (sub[0] & 0x1F000) >> 12;
        subevt[ind].elen = (sub[0] & 0x7FFE0000) >> 17;
        subevt[ind].fcode = (sub[0] & 0x80000000) >> 31;
        subevt[ind].time = ( (long long int)(sub[2] & 0xFFFF) << 32) + sub[1];
        subevt[ind].ctime = (sub[2] & 0x7FFF0000) >> 16;
        subevt[ind].ctimef = (sub[2] & 0x80000000) >> 31;
        subevt[ind].energy = (sub[3] & 0xFFFF);
        subevt[ind].trlen = (sub[3] & 0x7FFF0000) >> 16;
        subevt[ind].trwlen = subevt[ind].trlen / 2;
        subevt[ind].extra = (sub[3] & 0x80000000) >> 31;       
 
         
	float temp;
	temp=(float) subevt[ind].energy/2.;
	subevt[ind].energy=(int) temp;
	//if(subevt->trlen ==0){        
	//printf("chn: %d crn: %d sln: %d id: %d hlen: %d elen: %d energy: %d trlen: %d ribhold: %d ribh: %d fragsize: %d flagb: %d\n",subevt->chn,subevt->crn,subevt->sln,subevt->id,subevt->hlen,subevt->elen,subevt->energy,subevt->trlen,rib_sizeold,rib_size,frag_size,flagb);
	printf("chn: %d crn: %d sln: %d id: %d hlen: %d elen: %d ind: %d fragsize: %d temp: %d\n",subevt[ind].chn,subevt[ind].crn,subevt[ind].sln,subevt[ind].id, subevt[ind].hlen, subevt[ind].elen,ind,frag[ind].bsize,temporary);
       	//secraw(subevt);
	 //continue on if no trace, esum, or qsum
        if (subevt[ind].hlen==HEADER_LENGTH && subevt[ind].trwlen==0 ) {
           ind=ind+1;
		 continue;
        }
        //more data than just the header; read entire sub event
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        if (fread(sub, sizeof(int)*subevt[ind].elen, 1, fpr) != 1) break;
                              
        //trace
        k=0;
        for (i = subevt[ind].hlen; i < subevt[ind].elen; i++) {      
            subevt[ind].tr[i - subevt[ind].hlen + k] = sub[i] & 0x3FFF; // the upper 2 bits/16 bits are filled with 0s
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0xFFFF; //C. W  
           subevt[ind].tr[i - subevt[ind].hlen + k + 1] = (sub[i]>>16) & 0x3FFF;
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0xFFFF; //C. W
            k=k+1;
        } 


       	//printf("id: %d hlen: %d elen: %d\n",subevt->id,subevt->hlen,subevt->elen); 
	if(subevt[ind].trlen != 0){
	int l,m;
	for(m=0;m<8;m++){
	subevt[ind].qsum[m]=0;
	}
	for(l=0;l<31;l++){
	subevt[ind].qsum[0]=subevt[ind].qsum[0]+subevt[ind].tr[l];
	}
	for(l=31;l<60;l++){
	subevt[ind].qsum[1]=subevt[ind].qsum[1]+subevt[ind].tr[l];
	}
	for(l=60;l<75;l++){
	subevt[ind].qsum[2]=subevt[ind].qsum[2]+subevt[ind].tr[l];
	}
	for(l=75;l<95;l++){
	subevt[ind].qsum[3]=subevt[ind].qsum[3]+subevt[ind].tr[l];
	}
	for(l=95;l<105;l++){
	subevt[ind].qsum[4]=subevt[ind].qsum[4]+subevt[ind].tr[l];
	}
	for(l=105;l<160;l++){
	subevt[ind].qsum[5]=subevt[ind].qsum[5]+subevt[ind].tr[l];
	}
	for(l=160;l<175;l++){
	subevt[ind].qsum[6]=subevt[ind].qsum[6]+subevt[ind].tr[l];
	}
	for(l=175;l<200;l++){
	subevt[ind].qsum[7]=subevt[ind].qsum[7]+subevt[ind].tr[l];
	}
	}
     // if (subevt[sevtmult].id == 4 && subevt[sevtmult].fcode == 1) DB(subevt[sevtmult].tr);
            
        //continue if no esum or qsum   
        if (subevt[ind].hlen==HEADER_LENGTH) {
	   // pidevt(subevt);
		ind=ind+1;
            continue;
        }
        
        //esum
        if (subevt[ind].hlen==8 || subevt[ind].hlen==16) { 
            for (i=4; i < 8; i++) {
                subevt[ind].esum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[ind].hlen==12) { 
            for (i=4; i < 12; i++) {
                subevt[ind].qsum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[ind].hlen==16) { 
            for (i=8; i < 16; i++) {
                subevt[ind].qsum[i-8] = sub[i];
            }
        }    
     ind=ind+1;
           }//end loop over fragments for each ring
 	
	secraw(ind, subevt);
	
	}
}



