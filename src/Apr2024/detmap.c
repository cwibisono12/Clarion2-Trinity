//Detector Mapping:
//C. Wibisono
//Map the Raw subevent into corresponding types of Detector:
//Germanium mapping is refurbised from ORNL program scan.c
#include "global.h"
#include "mapping.h"
#include <math.h>
#include <stdbool.h>


void detmaps(int sevtmult, struct subevent *subevt, struct gdetector *ge, struct sidetector *si){
int i,j,l;     
 for (i=0; i<sevtmult; i++) { 
        
/* //07/31/'22 C.W

        for (j=0; j<sevtmult; j++) {
            if (i!=j)
            idid[subevt[i].id][subevt[j].id]++;
        }

*/
        
        //printf("i=%d, sevtmult=%d, subevt[i].id=%d, map2type[subevt[i].id]=%c, subevt[i].energy=%d, subevt[i].time=%lld\n", i, sevtmult, subevt[i].id, map2type[subevt[i].id], subevt[i].energy, subevt[i].time);
        //fflush(stdout);
          
        //if CFD is enabled, ctime will be non-zero
        //tempf = (float)subevt[i].ctime*10.0/32768.0;
        //tempf = (float)subevt[i].ctime*1.0/32768.0; //should make no difference ; need decimals or convert time to 1 ns bins
        //printf("%lld e=%d ", subevt[i].time, subevt[i].energy);
        //subevt[i].time = subevt[i].time + (long long int)tempf;  
        //printf("%lld and %f\n", subevt[i].time, tempf); 
  
               
        //Histogram calibrated tdif spectra (Do here to keep subevents within an event time ordered during event build above)
       /* 07/31/'22 //C.W
	 if (i==0) etime = subevt[i].time;
        tdif = abs(subevt[i].time - etime); 
        if (tdif >= 0 && tdif < 4096 && i!=0) tdif_cal[subevt[i].id][tdif]++;
      */
	  //if (tdif > 20 && tdif < 50) DB(subevt[i].tr);
        
        //tdif with respect to channel id 0
     /*  07/31/'22 C.2
	   if (i!=0 && subevt[0].id==0 && subevt[0].energy >= 10 && subevt[0].energy <= 8000 && subevt[i].energy >= 10  && subevt[i].energy <= 8000) {
            if (tdif >= 0 && tdif < 4096 ) tdif_cal0_ethresh[subevt[i].id][tdif]++;   
        }            
      */  

        /////////////////////
        // G = Ge Detector //
        /////////////////////              
        if ( map2type[subevt[i].id] == 'G' ) { //Keep G and ge or switch to C and clover/cl?
        
            if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_GE_XTL) { //Ge crystal 
               if ( ge[map2det[subevt[i].id]].xmult >= MAX_GE_XTL) {               
                    printf("SEVERE ERROR: Same Ge(xtl) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].xid[ge[map2det[subevt[i].id]].xmult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].xe[ge[map2det[subevt[i].id]].xmult] = subevt[i].energy; 
		ge[map2det[subevt[i].id]].xedopp[ge[map2det[subevt[i].id]].xmult] = subevt[i].energy/(1+map2beta[subevt[i].id]*cos(PIR*mapanglesi[subevt[i].id][0])); //C .W              
                ge[map2det[subevt[i].id]].xt[ge[map2det[subevt[i].id]].xmult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].xct[ge[map2det[subevt[i].id]].xmult] = subevt[i].ctime;  
                ge[map2det[subevt[i].id]].xpileup[ge[map2det[subevt[i].id]].xmult] = subevt[i].fcode;  
                ge[map2det[subevt[i].id]].xsubevtid[ge[map2det[subevt[i].id]].xmult] = i;  
                        
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][0] = mapangles[subevt[i].id][0];   
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][1] = mapanglesi[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][2] = mapangles1[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][3] = mapangles2[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][0] = mapangles[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][1] = mapanglesi[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][2] = mapangles1[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][3] = mapangles2[subevt[i].id][1]; 

                ge[map2det[subevt[i].id]].xmult++;   
            }                                                                
            if (map2deti[subevt[i].id] > MAX_GE_XTL && map2deti[subevt[i].id] <= MAX_GE_XTL + MAX_GE_SEG) { //Ge segment
               if ( ge[map2det[subevt[i].id]].smult >= MAX_GE_SEG ) {               
                    printf("SEVERE ERROR: Same Ge(seg) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].sid[ge[map2det[subevt[i].id]].smult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].se[ge[map2det[subevt[i].id]].smult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].st[ge[map2det[subevt[i].id]].smult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].sct[ge[map2det[subevt[i].id]].smult] = subevt[i].ctime;  
                ge[map2det[subevt[i].id]].spileup[ge[map2det[subevt[i].id]].smult] = subevt[i].fcode;  
                ge[map2det[subevt[i].id]].ssubevtid[ge[map2det[subevt[i].id]].smult] = i;  
                            
                ge[map2det[subevt[i].id]].smult++;               
            }            
            if (map2deti[subevt[i].id] > MAX_GE_XTL + MAX_GE_SEG) { //BGO 
               if ( ge[map2det[subevt[i].id]].bgomult >= MAX_GE_BGO ) {               
                    printf("SEVERE ERROR: Same Ge(bgo) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].bgoid[ge[map2det[subevt[i].id]].bgomult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].bgoe[ge[map2det[subevt[i].id]].bgomult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].bgot[ge[map2det[subevt[i].id]].bgomult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].bgoct[ge[map2det[subevt[i].id]].bgomult] = subevt[i].ctime;              
                ge[map2det[subevt[i].id]].bgopileup[ge[map2det[subevt[i].id]].bgomult] = subevt[i].fcode;
                ge[map2det[subevt[i].id]].bgosubevtid[ge[map2det[subevt[i].id]].bgomult] = i;  

                ge[map2det[subevt[i].id]].bgomult++;               
            } 
                 
                
        } //end G

        ///////////////////           //C. W
        // S=Si Detector //
	if ( map2type[subevt[i].id] == 'S' && subevt[i].trlen == 0 ){
	if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_SI_PAIR){
		            si[map2det[subevt[i].id]].siid[si[map2det[subevt[i].id]].simult] = map2deti[subevt[i].id];               
                si[map2det[subevt[i].id]].sie[si[map2det[subevt[i].id]].simult] = subevt[i].energy;             
                si[map2det[subevt[i].id]].sit[si[map2det[subevt[i].id]].simult] = subevt[i].time;                           
                si[map2det[subevt[i].id]].sict[si[map2det[subevt[i].id]].simult] = subevt[i].ctime;  
                si[map2det[subevt[i].id]].sipileup[si[map2det[subevt[i].id]].simult] = subevt[i].fcode;  
                si[map2det[subevt[i].id]].sisubevtid[si[map2det[subevt[i].id]].simult] = i;        
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][0];   
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][0];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][1]; 
		            si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] = subevt[i].qsum[0];   
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] = subevt[i].qsum[1];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] = subevt[i].qsum[2];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] = subevt[i].qsum[3];  
		            si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] = subevt[i].qsum[4];   
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] = subevt[i].qsum[5];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] = subevt[i].qsum[6];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] = subevt[i].qsum[7];              
                si[map2det[subevt[i].id]].simult++;   
                
		}

} //end Si        //C. W

// Another types of GAGG-mappings:
if ( map2type[subevt[i].id] == 'S' && subevt[i].trlen != 0 ){
	if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_SI_PAIR){
		            si[map2det[subevt[i].id]].siid[si[map2det[subevt[i].id]].simult] = map2deti[subevt[i].id];               
                si[map2det[subevt[i].id]].sie[si[map2det[subevt[i].id]].simult] = subevt[i].energy;             
                si[map2det[subevt[i].id]].sit[si[map2det[subevt[i].id]].simult] = subevt[i].time;                           
                si[map2det[subevt[i].id]].sict[si[map2det[subevt[i].id]].simult] = subevt[i].ctime;  
                si[map2det[subevt[i].id]].sipileup[si[map2det[subevt[i].id]].simult] = subevt[i].fcode;  
                si[map2det[subevt[i].id]].sisubevtid[si[map2det[subevt[i].id]].simult] = i;        
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][0];   
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][0];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][1]; 
//Reconstructing Trace to Qsum: //C.W
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] = 0;

for (l=0; l<31; l++){
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] + subevt[i].tr[l];
}   

for (l=31;l<60;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] + subevt[i].tr[l];                    
}

for (l=60;l<75;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] + subevt[i].tr[l];               
}

for (l=75;l<95;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] + subevt[i].tr[l];
}  

for (l=95;l<105;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] + subevt[i].tr[l];   
}

for (l=105;l<160;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] + subevt[i].tr[l];
}

for (l=160;l<175;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] + subevt[i].tr[l];               
}

for (l=175;l<200;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] + subevt[i].tr[l];              
}

si[map2det[subevt[i].id]].simult++;   
                
		}

} //end Si        //C. W


        ///////////////////  
        
      } // end i loop over sevtmult

}     

