//Unpack Pixie16 Digitizer and Build One Event 
//Based on work of J.M Allmond from scan.c
//Repackaged into a function by C. Wibisono

#include <stdio.h>
#include <stdlib.h>
#include "global.h"
#include "spec_dec.h"

long long int etime,tdif, idtime[MAX_ID]={0},temptime;
int toreader(unsigned int sub[], struct subevent *subevt, FILE *fpr){
      /////////////////////////////////
      // UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////
      int sevtmult;
      float tempf=0;
      int i=0,j=0,k=0;
	etime=-1; tdif=-1; sevtmult=0;  
      //memset(&subevt, 0, sizeof(subevt)); //not needed since everything is redefined (except maybe trace on pileup evts)
      while (1) { //get subevents and event build for one "event" 
        
       // memset(&subevt[sevtmult], 0, sizeof(subevt[sevtmult])); //not needed since everything is redefined (except maybe trace on pileup evts)
        
        //read 4-byte header
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        subevt[sevtmult].chn = sub[0] & 0xF;
        subevt[sevtmult].sln = (sub[0] & 0xF0) >> 4;
        subevt[sevtmult].crn = (sub[0] & 0xF00) >> 8;
        subevt[sevtmult].id = subevt[sevtmult].crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (subevt[sevtmult].sln-BOARD_START)*MAX_CHANNELS_PER_BOARD + subevt[sevtmult].chn;   
        subevt[sevtmult].hlen = (sub[0] & 0x1F000) >> 12;
        subevt[sevtmult].elen = (sub[0] & 0x7FFE0000) >> 17;
        subevt[sevtmult].fcode = (sub[0] & 0x80000000) >> 31;
        subevt[sevtmult].time = ( (long long int)(sub[2] & 0xFFFF) << 32) + sub[1];
        subevt[sevtmult].ctime = (sub[2] & 0x7FFF0000) >> 16;
        subevt[sevtmult].ctimef = (sub[2] & 0x80000000) >> 31;
        subevt[sevtmult].energy = (sub[3] & 0xFFFF);
        subevt[sevtmult].trlen = (sub[3] & 0x7FFF0000) >> 16;
        subevt[sevtmult].trwlen = subevt[sevtmult].trlen / 2;
        subevt[sevtmult].extra = (sub[3] & 0x80000000) >> 31;       
 
        //rebin raw trap energy from 32k to ....            
        tempf = (float)subevt[sevtmult].energy/RAWE_REBIN_FACTOR;// + RAND;       
        subevt[sevtmult].energy = (int)tempf;                        
 
        //check lengths (sometimes all of the bits for trace length are turned on ...)
       /* if (subevt[sevtmult].elen - subevt[sevtmult].hlen != subevt[sevtmult].trwlen) {
            printf("SEVERE ERROR: event, header, and trace length inconsistencies found\n");
            printf("event length = %d\n", subevt[sevtmult].elen);
            printf("header length = %d\n", subevt[sevtmult].hlen);
            printf("trace length = %d\n", subevt[sevtmult].trwlen);  
            printf("Extra = %d\n", subevt[sevtmult].extra); 
            printf("fcode = %d\n", subevt[sevtmult].fcode);              
            //sleep(1);          
            //return 0;
        } */ 
        
       
        //Set reference time for event building
        if (etime == -1) {
            etime = subevt[sevtmult].time;
            tdif = 0;
        }
        else {
            tdif = subevt[sevtmult].time - etime;
            if (tdif < 0) {
                printf("SEVERE ERROR: tdiff < 0, file must be time sorted\n");
                printf("etime = %lld, time = %lld, and tdif = %lld\n", etime, subevt[sevtmult].time, tdif);                
                return 0;   
            }    
        }    
      
        //Check for end of event, rewind, and break out of while loop
        if (tdif > EVENT_BUILD_TIME) {
            fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR); //fwrite/fread is buffered by system ; storing this in local buffer is no faster!
            break;           
        }    
               
        
        //time between sequential events for a single channel ; useful for determining optimal event building time
        temptime = (subevt[sevtmult].time - idtime[subevt[sevtmult].id])/100; //rebin to 1 micro-second
        if ( temptime >= 0 && temptime < 8192) {
            tdifid[subevt[sevtmult].id][temptime]++;
        }    
        idtime[subevt[sevtmult].id]=subevt[sevtmult].time; //store time for next subevent of channel    
    
        // total pileups
        if (subevt[sevtmult].fcode==1) {
            pileupcount++;
        }
        
        //Histogram raw spectra
        hit[0][subevt[sevtmult].id]++;
        if (subevt[sevtmult].fcode==1) 
        hit[1][subevt[sevtmult].id]++;     
  
        if (subevt[sevtmult].energy >= 0 && subevt[sevtmult].energy < 8192)
        e_raw[subevt[sevtmult].id][subevt[sevtmult].energy]++;
 
        if (subevt[sevtmult].time/1000000000 >= 0 && subevt[sevtmult].time/1000000000 < 8192) // rebin to 10 seconds
        tevt_raw[subevt[sevtmult].id][subevt[sevtmult].time/1000000000]++; // rebin to 10 seconds

        if (subevt[sevtmult].ctime >= 0 && subevt[sevtmult].ctime < 8192)
        tcfd_raw[subevt[sevtmult].id][subevt[sevtmult].ctime]++;

        if (tdif >= 0 && tdif < 4096 && sevtmult!=0)
        tdif_raw[subevt[sevtmult].id][tdif]++;


        //if CFD is enabled, ctime will be non-zero
        //tempf = (float)subevt[sevtmult].ctime*10.0/32768.0;
        //subevt[sevtmult].time = subevt[sevtmult].time + (long long int)tempf;    

        //Calibrate energy and time
        //Apr23 '24 prefactor 0.5 to accomodate 2024 setting
        tempf = ((float)0.5*subevt[sevtmult].energy*ecal[subevt[sevtmult].id][1] + ecal[subevt[sevtmult].id][0])/new_gain;// + RAND;       
        subevt[sevtmult].energy = (int)tempf; 
        //subevt[sevtmult].time += (long long int)tcal[subevt[sevtmult].id][0];       
	
        //Histogram calibrated spectra
        if (subevt[sevtmult].energy >= 0 && subevt[sevtmult].energy < 8192)
        e_cal[subevt[sevtmult].id][subevt[sevtmult].energy]++;        
        
        if (subevt[sevtmult].time/1000000000 >= 0 && subevt[sevtmult].time/1000000000 < 8192)        
        tevt_cal[subevt[sevtmult].id][subevt[sevtmult].time/1000000000]++;
         
        //continue on if no trace, esum, or qsum
        if (subevt[sevtmult].hlen==HEADER_LENGTH && subevt[sevtmult].trwlen==0 ) {
            sevtmult++;
            continue;
        }
        
        //more data than just the header; read entire sub event
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        if (fread(sub, sizeof(int)*subevt[sevtmult].elen, 1, fpr) != 1) break;
                              
        //trace
        k=0;
        for (i = subevt[sevtmult].hlen; i < subevt[sevtmult].elen; i++) {      
            subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0x3FFF; // the upper 2 bits/16 bits are filled with 0s
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0xFFFF; //C. W  
           subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0x3FFF;
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0xFFFF; //C. W
            k=k+1;
        } 
        
     // if (subevt[sevtmult].id == 4 && subevt[sevtmult].fcode == 1) DB(subevt[sevtmult].tr);
            
        //continue if no esum or qsum   
        if (subevt[sevtmult].hlen==HEADER_LENGTH) {
            sevtmult++;        
            continue;
        }
        
        //esum
        if (subevt[sevtmult].hlen==8 || subevt[sevtmult].hlen==16) { 
            for (i=4; i < 8; i++) {
                subevt[sevtmult].esum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[sevtmult].hlen==12) { 
            for (i=4; i < 12; i++) {
                subevt[sevtmult].qsum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[sevtmult].hlen==16) { 
            for (i=8; i < 16; i++) {
                subevt[sevtmult].qsum[i-8] = sub[i];
            }
        }    
        
        sevtmult++;
     
      } //end while loop for unpacking sub events and event building for one "event"
 
 /*  
    if (sevtmult==0) break; //end main WHILE LOOP when out of events 
      mult[0][sevtmult]++; //Histogram raw sub event multiplicity 
      sevtcount += sevtmult;
      evtcount++; //event-built number
 */

return sevtmult;
}



