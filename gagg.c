//Usage: GAGG Data Processing:
//C. Wibisono
#include "global.h"
#include "mapping.h"
#include <math.h>
#include <stdbool.h>

int gaggproc(struct sidetector *si){

int gaggvalid=0;
int sicount;
int i;
for (i=1;i<MAX_SI;i++){
	
  if (si[i].simult == 2 && (((si[i].siid[0] == 1 && si[i].siid[1] == 2)) || (si[i].siid[0] == 2 && si[i].siid[1] == 1)) && si[i].sipileup[0] == 0 && si[i].sipileup[1] == 0
&& (si[i].sit[0]-si[i].sit[1] + 2000 > 1990 && si[i].sit[0]-si[i].sit[1]  + 2000 < 2010 )
) 
          {
si[i].peak =  (((si[i].siqdc[0][3]-(20./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1])) + (si[i].siqdc[1][3] - (20./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/20.);
si[i].tail = (((si[i].siqdc[0][5]-(55./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1])) + (si[i].siqdc[1][5] - (55./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/20.);


si[i].traceint=(((si[i].siqdc[0][2]+si[i].siqdc[0][3]+si[i].siqdc[0][4]+si[i].siqdc[0][5]+si[i].siqdc[0][6]-(115./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1]))+(si[i].siqdc[1][2]+si[i].siqdc[1][3]+si[i].siqdc[1][4]+si[i].siqdc[1][5]+si[i].siqdc[1][6]-(115./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/30.);

si[i].tpratio=4000.*si[i].tail/si[i].peak;

si[i].id = i;
si[i].valid = 1;

if (si[i].id == gaid[i-1]){
si[i].energy=100.*(gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}


//SiPM angle assignment:
for (sicount=0;sicount<3;sicount++){
si[i].theta[sicount]=si[i].sitheta[0][sicount];
si[i].phi[sicount]=si[i].siphi[0][sicount];
//if (sicount == 0)
//printf("theta: %f phi: %f\n",si[i].theta[sicount],si[i].phi[sicount]);
}


//Selecting GAGG events:
bool ingate = false;
    if (i == 1) { ingate = bantest(si[1].tpratio,si[1].energy,polyX21,polyY21,GNMPTS[0]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 2) { ingate = bantest(si[2].tpratio,si[2].energy,polyX22,polyY22,GNMPTS[1]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 3) { ingate = bantest(si[3].tpratio,si[3].energy,polyX23,polyY23,GNMPTS[2]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 4) { ingate = bantest(si[4].tpratio,si[4].energy,polyX24,polyY24,GNMPTS[3]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 5) { ingate = bantest(si[5].tpratio,si[5].energy,polyX25,polyY25,GNMPTS[4]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 6) { ingate = bantest(si[6].tpratio,si[6].energy,polyX26,polyY26,GNMPTS[5]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 7) { ingate = bantest(si[7].tpratio,si[7].energy,polyX27,polyY27,GNMPTS[6]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 8) { ingate = bantest(si[8].tpratio,si[8].energy,polyX28,polyY28,GNMPTS[7]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 9) { ingate = bantest(si[9].tpratio,si[9].energy,polyX29,polyY29,GNMPTS[8]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 10) { ingate = bantest(si[10].tpratio,si[10].energy,polyX210,polyY210,GNMPTS[9]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 11) { ingate = bantest(si[11].tpratio,si[11].energy,polyX41,polyY41,GNMPTS[10]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 12) { ingate = bantest(si[12].tpratio,si[12].energy,polyX42,polyY42,GNMPTS[11]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 13) { ingate = bantest(si[13].tpratio,si[13].energy,polyX43,polyY43,GNMPTS[12]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 14) { ingate = bantest(si[14].tpratio,si[14].energy,polyX44,polyY44,GNMPTS[13]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 15) { ingate = bantest(si[15].tpratio,si[15].energy,polyX45,polyY45,GNMPTS[14]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 16) { ingate = bantest(si[16].tpratio,si[16].energy,polyX46,polyY46,GNMPTS[15]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 17) { ingate = bantest(si[17].tpratio,si[17].energy,polyX47,polyY47,GNMPTS[16]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 18) { ingate = bantest(si[18].tpratio,si[18].energy,polyX48,polyY48,GNMPTS[17]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 19) { ingate = bantest(si[19].tpratio,si[19].energy,polyX49,polyY49,GNMPTS[18]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 20) { ingate = bantest(si[20].tpratio,si[20].energy,polyX410,polyY410,GNMPTS[19]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 21) { ingate = bantest(si[21].tpratio,si[21].energy,polyX411,polyY411,GNMPTS[20]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 22) { ingate = bantest(si[22].tpratio,si[22].energy,polyX412,polyY412,GNMPTS[21]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 23) { ingate = bantest(si[23].tpratio,si[23].energy,polyX413,polyY413,GNMPTS[22]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 24) { ingate = bantest(si[24].tpratio,si[24].energy,polyX414,polyY414,GNMPTS[23]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 25) { ingate = bantest(si[25].tpratio,si[25].energy,polyX415,polyY415,GNMPTS[24]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 26) { ingate = bantest(si[26].tpratio,si[26].energy,polyX416,polyY416,GNMPTS[25]);
                  if (ingate ==1)gaggvalid++;}

//Transforming variable to another :
//Only select si event that matches condition above:


if (gaggvalid > 0 && ingate == 1){
  si[gaggvalid-1]=si[i];
}


//Time Correlation Between Germanium and GAGG:
/*
int gage;
int gev;
for (gage=1;gage<MAX_SI;gage++){
        if (si[gage].valid != 1) continue;
        si[gage].time=(si[gage].sit[0]+si[gage].sit[1])/2;
        for (gev=0;gev<gmult;gev++){
            int tdiffgage=(ge[gev].time-si[gage].time) + 2000;
            if ((ge[gev].energy > 0 && ge[gev].energy < 4096) && (tdiffgage > 0  && tdiffgage <4096)){
            gagg_ge_tdiff[ge[gev].energy][tdiffgage]++;}
                                   }     
            }

*/
 
//}//end loop over gevalid //C.W
		}// end Si condition for GAGG firing:
    
	}//end Si
return gaggvalid;
}

