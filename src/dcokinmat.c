//Usage: GG-Matrix for DCO (Spin-Changes)
//C. Wibisono

#include <stdio.h>
#include "global.h"
#include "dco.h"


void dcokinmat(int gmult, int gaggvalid, struct gdetector *ge){
int i,j,tdif;
//Pure gg-matrix without any gate:
      //gamma-gamma time and etrue 

//Pure gg-matrix with any gate:
      //gamma-gamma time and etrue 
   
 if (gaggvalid == 1 && gmult > 0){
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
                    tdif = ge[i].time - ge[j].time + 2000; 
                //time difference matrix between ge-ge within prompt region:    
/* 
               if (tdif >= 0 && tdif < 4096 && ge[i].etrue >= 0 && ge[i].etrue < 5000 && ge[j].etrue >= 0 && ge[j].etrue < 5000 ) {
                    if (ge[i].etrue < ge[j].etrue)
                    gg_tdif[ge[i].etrue][tdif]++;
                    else 
                    gg_tdif[ge[j].etrue][tdif]++;
                }   

*/
                //prompt gamma-gamma  
                //Only 90 deg dets populate x axis ;forward and backward 48.24 and 131.75 degrees populate y axis
               if (tdif >= 1940 && tdif <= 2060 ) {
                    if (ge[i].etrue >= 0 && ge[i].etrue < 5000 && ge[j].etrue >= 0 && ge[j].etrue < 5000
                        &&((ge[i].theta[0] >= 40 && ge[i].theta[0] <=50))
                        && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95) 
                                      ) 
                    gg_prompt_purefpkin[(int)ge[i].etrue][(int)ge[j].etrue]++;
		  
                    if (ge[i].etrue >= 0 && ge[i].etrue < 5000 && ge[j].etrue >= 0 && ge[j].etrue < 5000
                        &&( (ge[i].theta[0] >= 130 && ge[i].theta[0] <=140))
                        && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95)) 
                    gg_prompt_purebpkin[(int)ge[i].etrue][(int)ge[j].etrue]++;


                }  
              //non-prompt gamma-gamma (mult by 0.2555)
                if ( (tdif >= 1840 && tdif <= 1900) || (tdif >= 2100 && tdif <= 2160 ) ){ 

                    if (ge[i].etrue >= 0 && ge[i].etrue < 5000 && ge[j].etrue >= 0 && ge[j].etrue < 5000
                      &&((ge[i].theta[0] >=40 && ge[i].theta[0] <=50))
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95) 
                      ) 
                    gg_nprompt_purefpkin[(int)ge[i].etrue][(int)ge[j].etrue]++;
          
                    if (ge[i].etrue >= 0 && ge[i].etrue < 5000 && ge[j].etrue >= 0 && ge[j].etrue < 5000
                      &&( (ge[i].theta[0] >=130 && ge[i].theta[0] <=140))
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95) 
                      ) 
                    gg_nprompt_purebpkin[(int)ge[i].etrue][(int)ge[j].etrue]++;



                 }    
            }
        }
      }  

}

}
