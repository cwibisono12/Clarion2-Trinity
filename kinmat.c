//Kinematic Correction Function:
//C. Wibisono

#include "global.h"
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//Kinmeatic Correction:
void doppcorrect(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int option, float ml, float mbeam, float mh, float vbz, float betascale){
//Correcting Energy ://C.W
int gamrev;
if (gaggvalid == 1 && gmult > 0){
	//printf("%d %d\n",gaggvalid,gmult);
	for (gamrev=0;gamrev<gmult;gamrev++){
	//Start Doing Doppler Correction: //C.W Jun8:
	//Defining a Unit Vector:
	ge[gamrev].Dg[0]=sin((PIR)*ge[gamrev].theta[1])*cos((PIR)*ge[gamrev].phi[1]);
	ge[gamrev].Dg[1]=sin((PIR)*ge[gamrev].theta[1])*sin((PIR)*ge[gamrev].phi[1]);
	ge[gamrev].Dg[2]=cos((PIR)*ge[gamrev].theta[1]);

	//Deducing velocity of light particles:
	si[0].vl[0]=pow(2*(si[0].energy/100.)/ml,0.5)*sin((PIR)*si[0].theta[0])*cos((PIR)*si[0].phi[0]);
	si[0].vl[1]=pow(2*(si[0].energy/100.)/ml,0.5)*sin((PIR)*si[0].theta[0])*sin((PIR)*si[0].phi[0]);
	si[0].vl[2]=pow(2*(si[0].energy/100.)/ml,0.5)*cos((PIR)*si[0].theta[0]);


	//Deducing velocity of heavy particles:		
	si[0].vh[0]=-ml*si[0].vl[0]/mh;
	si[0].vh[1]=-ml*si[0].vl[1]/mh;
	si[0].vh[2]=(mbeam*vbz-ml*si[0].vl[2])/mh;

	si[0].betarest=pow(pow(si[0].vh[0],2)+pow(si[0].vh[1],2)+pow(si[0].vh[2],2),0.5);
	si[0].thetarest=100.0*acos((si[0].vh[2])/(si[0].betarest))*180.0/3.14159265;

	ge[gamrev].betavalue=(si[0].vh[0])*ge[gamrev].Dg[0]+(si[0].vh[1])*ge[gamrev].Dg[1]+(si[0].vh[2])*ge[gamrev].Dg[2];
	ge[gamrev].alpha=sin((PIR)*si[0].theta[0])*sin((PIR)*ge[gamrev].theta[1])*cos((PIR)*(si[0].phi[0]-ge[gamrev].phi[1]))+cos((PIR)*ge[gamrev].theta[1])*cos((PIR)*si[0].theta[0]);

//Non-relativistiv Doppler:
  if (option == 1) ge[gamrev].etrue=ge[gamrev].energy/(1+betascale*ge[gamrev].betavalue); //Doppler with kinmat
  else
	ge[gamrev].etrue=ge[gamrev].energy/(1+betascale*cos(PIR*ge[gamrev].theta[1])); //Doppler without Kinmat
//Relativistic Doppler Correction
	ge[gamrev].etrue2=ge[gamrev].energy*pow((1-pow(si[0].betarest,2)),0.5)/(1-(ge[gamrev].alpha)*(si[0].betarest));
//       printf("energync: %d,energyc: %f,betavalue1: %f, betavalue2: %f\n",ge[gamrev].energy,ge[gamrev].etrue,si[0].betarest,ge[gamrev].betavalue);

//Angle of Residual with respect to the beam axis:
						}			

//printf("%s , energy: %d, etrue: %d\n","doppcorrect",ge[0].energy,(int)ge[0].etrue);

					}
}
//==========================================



