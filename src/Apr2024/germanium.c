//Usage: Germanium Data Processing including Addback and Compton Supression
//C. Wibisono
//Refurbished from ONRL Program "scan.c"
#include "global.h"
#include "spec_dec.h"
#include <math.h>


int gamproc(struct gdetector *ge){
/////////////////////
// G = Ge Detector //
/////////////////////      
int i,j,k;
long long int tdif;
int max1,max2,maxid1,maxid2;
int gmult,gcount;

      gmult = 0;
      for (i=0; i<MAX_GE; i++) {
        
        max1 = -1; max2 = -1;
        maxid1 = -1; maxid2 = -1;
        
        ge[i].id = i;
        ge[i].xvalid =0;  //C.W
      //Addback and Compton Suppression 
	      for (j=0; j<ge[i].xmult; j++) {
            
            //compton suppression per crystal 
             if (GE_BGO_SUPPRESSION == TRUE) {
                for (k=0; k<ge[i].bgomult; k++) {
                    tdif = abs( ge[i].xt[j] - ge[i].bgot[k] );
                    if (tdif >= 0 && tdif <= 4096) ge_bgo_tdif[i][tdif]++;
                    if ( (tdif < 50 && ge[i].bgoe[k] > 10 && ge[i].bgoe[k] < 30000) || ge[i].bgopileup[k]==TRUE) { //need to fix bgo pileup with trace analysis 
                        ge[i].xsuppress[j] = TRUE; 
                        ge[i].suppress = TRUE;                        
                    }  
                }
            }      
            
            //addback     
            if (ge[i].xsuppress[j] == FALSE) {       
                if (ge[i].xpileup[j] == TRUE) {
                    ge[i].pileup = TRUE;     
                    continue;
                }              
                //xtl spectra
                if (ge[i].xe[j] > 0 && ge[i].xe[j] < 8192 && ge[i].id >= 1) ge_spe_xtl[(ge[i].id-1)*MAX_GE_XTL + ge[i].xid[j]][ge[i].xe[j]]++;                        
                tdif = abs( ge[i].xt[j] - ge[i].time );
      //          if (tdif >= 0 && tdif <= 4096 && ge[i].time != 0 ) ge_xtl_tdif[i][tdif]++;           
                                       
                if (ge[i].xe[j] > 50 && ge[i].xe[j] < 5000) {                  
    //                if (tdif >= 0 && tdif <= 4096 && ge[i].time != 0 ) ge_xtl_tdif_ethresh[i][tdif]++;           
                    if (tdif < 20 || ge[i].time == 0) {  
                        ge[i].energy = ge[i].energy + ge[i].xe[j];
			                  ge[i].xevalid[ge[i].xvalid]=ge[i].xe[j]; //C.W 
                        //ge[i].thetavalid[xvalid]=ge[i].xtheta[j][1]; //C.W
                        ge[i].xidvalid[ge[i].xvalid]=ge[i].xid[j]; //C.W
                        ge[i].xvalid++; //C.W
                        ge[i].edop = ge[i].edop + (int) ge[i].xedopp[j]; //C.W
                        if (max1 < ge[i].xe[j]) {
                            max1 = ge[i].xe[j];
                            maxid1 = j;
                            ge[i].time = ge[i].xt[j];
                            ge[i].ctime = ge[i].xct[j];
                        } 
                    }
                    else {
                        ge[i].nonprompt = TRUE; // the first time will become the adopted value / event   
                    }        
                } 
            }//end for ge.xsuppress  --C. W  
           
	}//enf for germanium crystal multiplicity --C. W
        
        if (max1 == -1) continue;
        
        //Segmentation Position and Compton Suppression
        for (j=0; j<ge[i].smult; j++) {

            //compton suppression per segment 
            if (GE_BGO_SUPPRESSION == TRUE) {            
                for (k=0; k<ge[i].bgomult; k++) {
                    tdif = abs( ge[i].st[j] - ge[i].bgot[k] );
                    if (tdif < 50 && ge[i].bgoe[k] > 10) {
                        ge[i].ssuppress[j] = TRUE;    
                    }
                }
            }      

            //segment
            if (ge[i].ssuppress[j] == FALSE && ge[i].se[j] > 0 && ge[i].se[j] < 10000) {                    
                if (max2 < ge[i].se[j]) {
                    max2 = ge[i].se[j];
                    maxid2 = j;
                }                                           
            }  

        }    
        
        //Angle assignments
        ge[i].theta[0] = ge[i].xtheta[maxid1][0];                       //detector center
        ge[i].phi[0] = ge[i].xphi[maxid1][0];
        ge[i].theta[1] = ge[i].xtheta[maxid1][1];                       //crystal center
        ge[i].phi[1] = ge[i].xphi[maxid1][1];       
        if (ge[i].sid[maxid2] == 6) {                                   // side channel C  
            ge[i].theta[2] = ge[i].xtheta[maxid1][2];
            ge[i].phi[2] = ge[i].xphi[maxid1][2];  
        }
        else if (ge[i].sid[maxid2] == 5 || ge[i].sid[maxid2] == 7) {    // side channel L/R
            ge[i].theta[2] = ge[i].xtheta[maxid1][3];
            ge[i].phi[2] = ge[i].xphi[maxid1][3];             
        }  
        else {      
            ge[i].theta[2] = ge[i].xtheta[maxid1][1];                   // side channel failure --> crystal center
            ge[i].phi[2] = ge[i].xphi[maxid1][1];  
        }        
        
        //clean addback
//        if (ge[i].suppress == FALSE && ge[i].pileup == FALSE && ge[i].nonprompt == FALSE) {
        if (ge[i].suppress == FALSE && ge[i].nonprompt == FALSE) {

            ge[i].clean=TRUE; //maybe clean should not include suppress == FALSE?
        }    
        

        //Ge spectra
        if (ge[i].energy > 0 && ge[i].energy < 8192){ 
	    ge_spe[ge[i].id][ge[i].energy]++;
//	    ge_spe_dopp[ge[i].id][ge[i].edop]++; //C.W
	}         
  //      if (ge[i].energy > 0 && ge[i].energy < 4096 && ge[i].clean == TRUE) ge_spe_clean[ge[i].id][ge[i].energy]++;
        
        //copy data and increment counters
        ge[gmult]=ge[i];
        gmult++;        
        gcount++;
        
      } //end G
return gmult;
}






