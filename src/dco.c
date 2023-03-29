//Usage: GG-Matrix for DCO (Spin-Changes)
//C. Wibisono

#include <stdio.h>
#include "global.h"
#include "dco.h"


void dco(int gmult, int gaggvalid, struct gdetector *ge){
int i,j,tdif;
//Pure gg-matrix without any gate:
      //gamma-gamma time and energy 
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
                    tdif = ge[i].time - ge[j].time + 2000; 
 
               //time difference matrix between ge-ge within prompt region:    
/*
                if (tdif >= 0 && tdif < 4096 && ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000 ) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif[ge[i].energy][tdif]++;
                    else 
                    gg_tdif[ge[j].energy][tdif]++;
                }   

*/
                //prompt gamma-gamma  
                //Only 90 deg dets populate x axis ;forward and backward 48.24 and 131.75 degrees populate y axis
               if (tdif >= 1940 && tdif <= 2060 ) {
                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                        &&((ge[i].theta[0] >= 40 && ge[i].theta[0] <=50))
                        && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95) 
                                      ) 
                    gg_prompt_puref[ge[i].energy][ge[j].energy]++;
		  
                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                        &&( (ge[i].theta[0] >= 130 && ge[i].theta[0] <=140))
                        && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95)) 
                    gg_prompt_pureb[ge[i].energy][ge[j].energy]++;

/*
             if(ge[i].etrue2 >=0 && ge[i].etrue2 < 5000 && ge[j].energy >=0 && ge[j].energy < 5000
                      && ((ge[i].theta[0] >= 40 && ge[i].theta[0] <= 50) || (ge[i].theta[0] >= 130 && ge[i].theta[0] <= 140))
                      && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95)
){
               gg_prompt_nddopp_pure[(int)ge[i].etrue2][ge[j].energy]++;} //C.W
              if (ge[i].etrue2 >=0 && ge[i].etrue2 < 5000 && ge[j].etrue2 >=0 && ge[j].etrue2 <5000
                       && ((ge[i].theta[0] >= 40 && ge[i].theta[0] <=50) || (ge[i].theta[0] >=130 && ge[i].theta[0] <=140) )
                       && (ge[j].theta[0] >= 85 && ge[j].theta[0] <= 95)
                              ){
               gg_prompt_ddopp_pure[(int)ge[i].etrue2][(int)ge[j].etrue2]++;} //C.W                   
            */
                }  
              //non-prompt gamma-gamma (mult by 0.2555)
                if ( (tdif >= 1840 && tdif <= 1900) || (tdif >= 2100 && tdif <= 2160 ) ){ 

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                      &&((ge[i].theta[0] >=40 && ge[i].theta[0] <=50))
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95) 
                      ) 
                    gg_nprompt_puref[ge[i].energy][ge[j].energy]++;
          
                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                      &&( (ge[i].theta[0] >=130 && ge[i].theta[0] <=140))
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95) 
                      ) 
                    gg_nprompt_pureb[ge[i].energy][ge[j].energy]++;


	/*
              if(ge[i].etrue2 >=0 && ge[i].etrue2 <5000 && ge[j].energy >=0 && ge[j].energy <5000
                      && ((ge[i].theta[0] >=40 && ge[i].theta[0] <= 50) || (ge[j].theta[0] >= 130 && ge[j].theta[0] <=140))
                      && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95)
        ){
            		    gg_nprompt_nddopp_pure[(int)ge[i].etrue2][ge[j].energy]++;} //C.W 
                if(ge[i].etrue2 >=0 && ge[i].etrue2 <5000 && ge[j].etrue2 >=0 && ge[j].etrue2 <5000
                      && ((ge[i].theta[0]>= 40 && ge[i].theta[0] <=50) || (ge[i].theta[0] >=130 && ge[i].theta[0] <= 140)) 
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95)
){ 
                   gg_nprompt_ddopp_pure[(int)ge[i].etrue2][(int)ge[j].etrue2]++;} //C.W              
               
            */
                 }    
            }
        }
      }  

//Pure gg-matrix with any gate:
      //gamma-gamma time and energy 
   
 if (gaggvalid == 1 && gmult > 0){
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
                    tdif = ge[i].time - ge[j].time + 2000; 
                //time difference matrix between ge-ge within prompt region:    
/* 
               if (tdif >= 0 && tdif < 4096 && ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000 ) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif[ge[i].energy][tdif]++;
                    else 
                    gg_tdif[ge[j].energy][tdif]++;
                }   

*/
                //prompt gamma-gamma  
                //Only 90 deg dets populate x axis ;forward and backward 48.24 and 131.75 degrees populate y axis
               if (tdif >= 1940 && tdif <= 2060 ) {
                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                        &&((ge[i].theta[0] >= 40 && ge[i].theta[0] <=50))
                        && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95) 
                                      ) 
                    gg_prompt_purefp[ge[i].energy][ge[j].energy]++;
		  
                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                        &&( (ge[i].theta[0] >= 130 && ge[i].theta[0] <=140))
                        && (ge[j].theta[0] >= 85 && ge[j].theta[0] <=95)) 
                    gg_prompt_purebp[ge[i].energy][ge[j].energy]++;


                }  
              //non-prompt gamma-gamma (mult by 0.2555)
                if ( (tdif >= 1840 && tdif <= 1900) || (tdif >= 2100 && tdif <= 2160 ) ){ 

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                      &&((ge[i].theta[0] >=40 && ge[i].theta[0] <=50))
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95) 
                      ) 
                    gg_nprompt_purefp[ge[i].energy][ge[j].energy]++;
          
                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000
                      &&( (ge[i].theta[0] >=130 && ge[i].theta[0] <=140))
                      && (ge[j].theta[0] >=85 && ge[j].theta[0] <=95) 
                      ) 
                    gg_nprompt_purebp[ge[i].energy][ge[j].energy]++;



                 }    
            }
        }
      }  

}

}
