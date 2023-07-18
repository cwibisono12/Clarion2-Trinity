//Polarization Function:
//C. Wibisono

#include "global.h"
#include "pol.h"
#include <math.h>
#include <stdbool.h>



//Polarization Analysis:
void pol(int gaggvalid, int gmult, struct gdetector *ge){

//UnGated:
if(gmult>0){
int a,b;
//for(b=0;b<gmult;b++){
	for(a=0;a<gmult;a++){
//	if(b!=a){	 
       if(ge[a].xvalid == 2 && (ge[a].id == 3 || ge[a].id == 4  || ge[a].id == 9)){
                 //Parallel Spectra:
                  if((ge[a].xidvalid[0] == 1 && ge[a].xidvalid[1] == 4) || (ge[a].xidvalid[0] == 4 && ge[a].xidvalid[1] == 1)
                    || (ge[a].xidvalid[0] == 2 && ge[a].xidvalid[1] == 3) || (ge[a].xidvalid[0] == 3 && ge[a].xidvalid[1] == 2)){
                      ge[a].epol=ge[a].xevalid[0]+ge[a].xevalid[1];
                     // printf("%f %d %d %d  %d\n",ge[a].theta[0],ge[a].xidvalid[0],ge[a].xidvalid[1],ge[a].xvalid,ge[a].epol);
                  //Histogram 1D spectra:
                   if((ge[a].epol > 0 && ge[a].epol < 5000))
                      pol_para[0][ge[a].epol]++;
                                                                   }
                
                  //Perpendicular Spectra:
                  if((ge[a].xidvalid[0] == 1 && ge[a].xidvalid[1] == 2) || (ge[a].xidvalid[0] == 2 && ge[a].xidvalid[1] == 1) 
                    || (ge[a].xidvalid[0] == 3 && ge[a].xidvalid[1] == 4) || (ge[a].xidvalid[0] == 4 && ge[a].xidvalid[1] == 3)){
                      ge[a].epol=ge[a].xevalid[0]+ge[a].xevalid[1];
                  //Histogram 1D spectra:
                     // printf("%f %d %d %d  %d\n",ge[a].theta[0],ge[a].xidvalid[0],ge[a].xidvalid[1],ge[a].xvalid,ge[a].epol);
                   if((ge[a].epol > 0 && ge[a].epol <5000))
                      pol_perp[0][ge[a].epol]++;
                                                                  }
                                            }

              		}      
	//	}
	//}
}


//Particle Gated:
if(gaggvalid == 1 && gmult >0){
int a,b;
//for(b=0;b<gmult;b++){
	for(a=0;a<gmult;a++){
//	if(b!=a){	 
       if(ge[a].xvalid == 2 && (ge[a].id == 3 || ge[a].id == 4  || ge[a].id == 9)){
                 //Parallel Spectra:
                  if((ge[a].xidvalid[0] == 1 && ge[a].xidvalid[1] == 4) || (ge[a].xidvalid[0] == 4 && ge[a].xidvalid[1] == 1)
                    || (ge[a].xidvalid[0] == 2 && ge[a].xidvalid[1] == 3) || (ge[a].xidvalid[0] == 3 && ge[a].xidvalid[1] == 2)){
                      ge[a].epol=ge[a].xevalid[0]+ge[a].xevalid[1];
                     // printf("%f %d %d %d  %d\n",ge[a].theta[0],ge[a].xidvalid[0],ge[a].xidvalid[1],ge[a].xvalid,ge[a].epol);
                  //Histogram 1D spectra:
                   if((ge[a].epol > 0 && ge[a].epol < 5000))
                      pol_para_p[0][ge[a].epol]++;
                                                                   }
                
                  //Perpendicular Spectra:
                  if((ge[a].xidvalid[0] == 1 && ge[a].xidvalid[1] == 2) || (ge[a].xidvalid[0] == 2 && ge[a].xidvalid[1] == 1) 
                    || (ge[a].xidvalid[0] == 3 && ge[a].xidvalid[1] == 4) || (ge[a].xidvalid[0] == 4 && ge[a].xidvalid[1] == 3)){
                      ge[a].epol=ge[a].xevalid[0]+ge[a].xevalid[1];
                  //Histogram 1D spectra:
                     // printf("%f %d %d %d  %d\n",ge[a].theta[0],ge[a].xidvalid[0],ge[a].xidvalid[1],ge[a].xvalid,ge[a].epol);
                   if((ge[a].epol > 0 && ge[a].epol <5000))
                      pol_perp_p[0][ge[a].epol]++;
                                                                  }
                                            }

              		}      
	//	}
	//}
}






}


