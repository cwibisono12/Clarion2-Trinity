//Create gamma-gamma coincidences with particle gated
//C. Wibisono

//Only gg that come in coincidences with proton allowed

#include "global.h"
#include "mapping.h"
#include "ggmat.h"
#include <math.h>
#include <stdio.h>

void ggmat(int gaggvalid, int gmult, struct gdetector *ge){
if(gaggvalid == 1 && gmult > 0){
long long int tdif1,tdif2;
int i,j;
      //gamma-gamma time and energy 
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
		/*    
                if (ge[i].validp == 1 && ge[j].validp == 1)
                    tdif1 = ge[i].time - ge[j].time + 2000;
                if (ge[i].validnp == 1 && ge[j].validnp == 1)
                    tdif2 = ge[i].time - ge[j].time + 2000; 
                //time difference matrix between ge-ge within prompt region:    
                if (tdif1 >= 0 && tdif1 < 4096 && ge[i].energy >= 0 && ge[i].energy < 4096 && ge[j].energy >= 0 && ge[j].energy < 4096 && (ge[i].validp == 1) && ge[j].validp == 1) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif1[ge[i].energy][tdif1]++;
                    else 
                    gg_tdif1[ge[j].energy][tdif1]++;
                }   
               //time difference matrix between ge-ge within non-prompt region:
                if (tdif2 >=0 && tdif2 < 4096 && ge[i].energy >=0 && ge[j].energy < 4096 && ge[j].energy >=0 && ge[j].energy < 4096 && ge[i].validnp == 1 && ge[j].validnp == 1) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif2[ge[i].energy][tdif2]++;
                    else
                    gg_tdif2[ge[j].energy][tdif2]++;
                }

		*/
                //prompt gamma-gamma
               // if (tdif1 >= 1967 && tdif1 <= 2032 && ge[i].validp == 1 && ge[j].validp == 1) {
                if (ge[i].validp == 1 && ge[j].validp == 1) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_prompt[ge[i].energy][ge[j].energy]++;
		          if(ge[i].etrue >=0 && ge[i].etrue < 5000 && ge[j].energy >=0 && ge[j].energy < 5000){
               gg_prompt_nddopp[(int)ge[i].etrue][ge[j].energy]++;} //C.W
              if (ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){
               gg_prompt_ddopp[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W
                }        
                
                //prompt-non prompt:
                
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validp == 1 && ge[j].validnp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt1[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp1[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
              if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp1[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validnp == 1 && ge[j].validp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt2[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp2[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp2[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                
                
                
                //non-prompt and prompt:
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validnp == 1 && ge[j].validnp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                           
            }
        }
      }  
}//end gagg_valid

if(gmult > 0){
long long int tdif1;
int i,j;
      //gamma-gamma time and energy 
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
               // if (ge[i].validp == 1 && ge[j].validp == 1)
                    tdif1 = ge[i].time - ge[j].time + 2000;
               // if (ge[i].validnp == 1 && ge[j].validnp == 1)
                //    tdif2 = ge[i].time - ge[j].time + 2000; 
                //time difference matrix between ge-ge within prompt region:    
                if (tdif1 >= 0 && tdif1 < 4096 && ge[i].energy >= 0 && ge[i].energy < 4096 && ge[j].energy >= 0 && ge[j].energy < 4096 && (ge[i].validp == 1) && ge[j].validp == 1) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif1[ge[i].energy][tdif1]++;
                    else 
                    gg_tdif1[ge[j].energy][tdif1]++;
                }   
               //time difference matrix between ge-ge within non-prompt region:
               /*
	       	if (tdif2 >=0 && tdif2 < 4096 && ge[i].energy >=0 && ge[j].energy < 4096 && ge[j].energy >=0 && ge[j].energy < 4096 && ge[i].validnp == 1 && ge[j].validnp == 1) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif2[ge[i].energy][tdif2]++;
                    else
                    gg_tdif2[ge[j].energy][tdif2]++;
                }
		*/

                //prompt gamma-gamma
                if (tdif1 >= 1967 && tdif1 <= 2032) {
               // if (ge[i].validp == 1 && ge[j].validp == 1) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_prompt_ng[ge[i].energy][ge[j].energy]++;
		          if(ge[i].etrue >=0 && ge[i].etrue < 5000 && ge[j].energy >=0 && ge[j].energy < 5000){
               gg_prompt_nddopp_ng[(int)ge[i].etrue][ge[j].energy]++;} //C.W
              if (ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){
               gg_prompt_ddopp_ng[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W
                }        
                
		/*
                //prompt-non prompt:
                
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validp == 1 && ge[j].validnp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt1[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp1[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
              if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp1[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validnp == 1 && ge[j].validp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt2[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp2[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp2[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                
                
                */
                //non-prompt and prompt:
                 //non-prompt gamma-gamma (mult by 0.2555)
                //if ( (tdif2 >= 1963 && tdif2 <= 2040) ) {
                if ( (tdif1 >= 1840 && tdif1 <= 1872) || (tdif1 >= 2130 && tdif1 <= 2162) ) {
               // if ( ge[i].validnp == 1 && ge[j].validnp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt_ng[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            	   gg_nprompt_nddopp_ng[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp_ng[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                           
            }
        }
      }  
}//end gagg_valid

}
