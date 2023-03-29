//Analysis of Gamma+GAGG:
//C. Wibisono
//Generate 1D-Histogram:
#include "global.h"
#include "spec_dec.h"
#include <math.h>
#include <stdbool.h>

void gammagagg1Dspectra(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge){
int gem;
int gaggi;
//1 Dimensional Histogram //C. W
if (gaggvalid == 1  && gmult > 0){
	//printf("energy: %d, etrue: %d\n",ge[0].energy, (int)ge[0].etrue);
        for (gem=0;gem<gmult;gem++){
           for (gaggi=0;gaggi<gaggvalid;gaggi++){
//put condition below if we want to have only one peak for the timediff spectra     
// if (si[gaggi].energy > 200 ){  
          si[gaggi].time=(si[gaggi].sit[0]+si[gaggi].sit[1])/2.;
         int tdiffgemgam=si[gaggi].time-ge[gem].time + 2000;
         if ((tdiffgemgam > 0 && tdiffgemgam < 4096)&&(si[gaggi].energy > 0 && si[gaggi].energy < 4096)){    
   	if (si[gaggi].id >= 1 && si[gaggi].id <=10)	
	    gagg_gamma_tdiff_ring2[(int)si[gaggi].energy][tdiffgemgam]++; 
	if (si[gaggi].id >=11 && si[gaggi].id <=26)
	    gagg_gamma_tdiff_ring4[(int)si[gaggi].energy][tdiffgemgam]++;
		}
         if ((tdiffgemgam > 0 && tdiffgemgam < 4096) && (ge[gem].energy > 0 && ge[gem].energy < 5000))
            gagg_gamma_tdiff2[ge[gem].energy][tdiffgemgam]++;
              
//Prompt gamma:
         // if ((ge[gem].energy > 0 && ge[gem].energy < 4096) && (tdiffgemgam > 1963 && tdiffgemgam < 2002)){
         if ((ge[gem].energy > 0 && ge[gem].energy < 5000) && (tdiffgemgam >=1960 && tdiffgemgam <=2002)){
	  gated_gamma[0][ge[gem].energy]++; 
          gated_gammatrue[0][(int)ge[gem].etrue]++;
	  ge[gem].validp = 1;
}
//Non-Prompt gamma:
       //   if ((ge[gem].energy > 0 && ge[gem].energy < 4096) &&( (tdiffgemgam > 2010 && tdiffgemgam < 2049))){
         if ((ge[gem].energy > 0 && ge[gem].energy < 5000) &&( (tdiffgemgam >= 1891 && tdiffgemgam <=1912) || (tdiffgemgam >= 2050 && tdiffgemgam <= 2071) )){
	  gated_gamma_np[0][ge[gem].energy]++;
	  gated_gamma_nptrue[0][(int)ge[gem].etrue]++;   
          ge[gem].validnp = 1;
            }                           
                            //  }
                    }
            }
}
}
void parthit(int gaggvalid, int gmult, struct gdetector *ge, struct sidetector *si){
//Histogram for Number of Hits:
if (gaggvalid >0 && gmult > 0){
if (gaggvalid == 1) {gagghist[0][1000]++;
if (si[0].id >= 1 && si[0].id <= 10){
	theta_res2[si[0].id][(int)si[0].thetarest]++;
	e_lit2[si[0].id][(int)si[0].energy]++;
//	printf("%d %d %d\n",gaggvalid,gmult,(int)si[0].thetarest);
}

//if (si[0].id >= 11) printf("%d\n",si[0].id);
if (si[0].id >= 11 && si[0].id <= 26){
	theta_res4[si[0].id][(int)si[0].thetarest]++;
	e_lit4[si[0].id][(int)si[0].energy]++;
//	printf("%d %d\n",gaggvalid,(int)si[0].thetarest);
}
}
if (gaggvalid == 2) gagghist[0][2000]++;
if (gaggvalid > 2) gagghist[0][3000]++;
}
}



