//Usage: GAGG Data Processing:
//C. Wibisono
#include "global.h"
#include "mapping.h"
#include <math.h>
#include <stdbool.h>

int gaggproc(struct sidetector *si, int parttype){
//printf("%d\n",parttype);
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


si[i].traceint=(((si[i].siqdc[0][2]+si[i].siqdc[0][3]+si[i].siqdc[0][4]+si[i].siqdc[0][5]+si[i].siqdc[0][6]-(115./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1]))+(si[i].siqdc[1][2]+si[i].siqdc[1][3]+si[i].siqdc[1][4]+si[i].siqdc[1][5]+si[i].siqdc[1][6]-(115./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/50.);

si[i].tpratio=4000.*si[i].tail/si[i].peak;

si[i].id = i;
si[i].valid = 1;

if (si[i].id == gaid[i-1]){
if (parttype==1){
si[i].energy=100.*(gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}
if (parttype==2){
si[i].energy=si[i].energy=100.*(gaggquad[i-1]*pow(si[i].traceint,2.0)+gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}
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
    if (i == 1) { ingate = bantest(si[1].tpratio,si[1].traceint,polyX11,polyY11,GNMPTS[0]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 2) { ingate = bantest(si[2].tpratio,si[2].traceint,polyX12,polyY12,GNMPTS[1]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 3) { ingate = bantest(si[3].tpratio,si[3].traceint,polyX13,polyY13,GNMPTS[2]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 4) { ingate = bantest(si[4].tpratio,si[4].traceint,polyX14,polyY14,GNMPTS[3]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 5) { ingate = bantest(si[5].tpratio,si[5].traceint,polyX15,polyY15,GNMPTS[4]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 6) { ingate = bantest(si[6].tpratio,si[6].traceint,polyX16,polyY16,GNMPTS[5]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 7) { ingate = bantest(si[7].tpratio,si[7].traceint,polyX17,polyY17,GNMPTS[6]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 8) { ingate = bantest(si[8].tpratio,si[8].traceint,polyX18,polyY18,GNMPTS[7]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 9) { ingate = bantest(si[9].tpratio,si[9].traceint,polyX21,polyY21,GNMPTS[8]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 10) { ingate = bantest(si[10].tpratio,si[10].traceint,polyX22,polyY22,GNMPTS[9]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 11) { ingate = bantest(si[11].tpratio,si[11].traceint,polyX23,polyY23,GNMPTS[10]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 12) { ingate = bantest(si[12].tpratio,si[12].traceint,polyX24,polyY24,GNMPTS[11]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 13) { ingate = bantest(si[13].tpratio,si[13].traceint,polyX25,polyY25,GNMPTS[12]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 14) { ingate = bantest(si[14].tpratio,si[14].traceint,polyX26,polyY26,GNMPTS[13]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 15) { ingate = bantest(si[15].tpratio,si[15].traceint,polyX27,polyY27,GNMPTS[14]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 16) { ingate = bantest(si[16].tpratio,si[16].traceint,polyX28,polyY28,GNMPTS[15]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 17) { ingate = bantest(si[17].tpratio,si[17].traceint,polyX29,polyY29,GNMPTS[16]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 18) { ingate = bantest(si[18].tpratio,si[18].traceint,polyX210,polyY210,GNMPTS[17]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 19) { ingate = bantest(si[19].tpratio,si[19].traceint,polyX31,polyY31,GNMPTS[18]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 20) { ingate = bantest(si[20].tpratio,si[20].traceint,polyX32,polyY32,GNMPTS[19]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 21) { ingate = bantest(si[21].tpratio,si[21].traceint,polyX33,polyY33,GNMPTS[20]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 22) { ingate = bantest(si[22].tpratio,si[22].traceint,polyX34,polyY34,GNMPTS[21]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 23) { ingate = bantest(si[23].tpratio,si[23].traceint,polyX35,polyY35,GNMPTS[22]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 24) { ingate = bantest(si[24].tpratio,si[24].traceint,polyX36,polyY36,GNMPTS[23]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 25) { ingate = bantest(si[25].tpratio,si[25].traceint,polyX37,polyY37,GNMPTS[24]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 26) { ingate = bantest(si[26].tpratio,si[26].traceint,polyX38,polyY38,GNMPTS[25]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 27) { ingate = bantest(si[27].tpratio,si[27].traceint,polyX39,polyY39,GNMPTS[26]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 28) { ingate = bantest(si[28].tpratio,si[28].traceint,polyX310,polyY310,GNMPTS[27]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 29) { ingate = bantest(si[29].tpratio,si[29].traceint,polyX311,polyY311,GNMPTS[28]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 30) { ingate = bantest(si[30].tpratio,si[30].traceint,polyX312,polyY312,GNMPTS[29]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 31) { ingate = bantest(si[31].tpratio,si[31].traceint,polyX313,polyY313,GNMPTS[30]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 32) { ingate = bantest(si[32].tpratio,si[32].traceint,polyX314,polyY314,GNMPTS[31]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 33) { ingate = bantest(si[33].tpratio,si[33].traceint,polyX41,polyY41,GNMPTS[32]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 34) { ingate = bantest(si[34].tpratio,si[34].traceint,polyX42,polyY42,GNMPTS[33]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 35) { ingate = bantest(si[35].tpratio,si[35].traceint,polyX43,polyY43,GNMPTS[34]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 36) { ingate = bantest(si[36].tpratio,si[36].traceint,polyX44,polyY44,GNMPTS[35]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 37) { ingate = bantest(si[37].tpratio,si[37].traceint,polyX45,polyY45,GNMPTS[36]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 38) { ingate = bantest(si[38].tpratio,si[38].traceint,polyX46,polyY46,GNMPTS[37]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 39) { ingate = bantest(si[39].tpratio,si[39].traceint,polyX47,polyY47,GNMPTS[38]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 40) { ingate = bantest(si[40].tpratio,si[40].traceint,polyX48,polyY48,GNMPTS[39]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 41) { ingate = bantest(si[41].tpratio,si[41].traceint,polyX49,polyY49,GNMPTS[40]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 42) { ingate = bantest(si[42].tpratio,si[42].traceint,polyX410,polyY410,GNMPTS[41]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 43) { ingate = bantest(si[43].tpratio,si[43].traceint,polyX411,polyY411,GNMPTS[42]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 44) { ingate = bantest(si[44].tpratio,si[44].traceint,polyX412,polyY412,GNMPTS[43]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 45) { ingate = bantest(si[45].tpratio,si[45].traceint,polyX413,polyY413,GNMPTS[44]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 46) { ingate = bantest(si[46].tpratio,si[46].traceint,polyX414,polyY414,GNMPTS[45]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 47) { ingate = bantest(si[47].tpratio,si[47].traceint,polyX415,polyY415,GNMPTS[46]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 48) { ingate = bantest(si[48].tpratio,si[48].traceint,polyX416,polyY416,GNMPTS[47]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 49) { ingate = bantest(si[49].tpratio,si[49].traceint,polyX51,polyY51,GNMPTS[48]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 50) { ingate = bantest(si[50].tpratio,si[50].traceint,polyX52,polyY52,GNMPTS[49]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 51) { ingate = bantest(si[51].tpratio,si[51].traceint,polyX53,polyY53,GNMPTS[50]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 52) { ingate = bantest(si[52].tpratio,si[52].traceint,polyX54,polyY54,GNMPTS[51]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 53) { ingate = bantest(si[53].tpratio,si[53].traceint,polyX55,polyY55,GNMPTS[52]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 54) { ingate = bantest(si[54].tpratio,si[54].traceint,polyX56,polyY56,GNMPTS[53]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 55) { ingate = bantest(si[55].tpratio,si[55].traceint,polyX57,polyY57,GNMPTS[54]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 56) { ingate = bantest(si[56].tpratio,si[56].traceint,polyX58,polyY58,GNMPTS[55]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 57) { ingate = bantest(si[57].tpratio,si[57].traceint,polyX59,polyY59,GNMPTS[56]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 58) { ingate = bantest(si[58].tpratio,si[58].traceint,polyX510,polyY510,GNMPTS[57]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 59) { ingate = bantest(si[59].tpratio,si[59].traceint,polyX511,polyY511,GNMPTS[58]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 60) { ingate = bantest(si[60].tpratio,si[60].traceint,polyX512,polyY512,GNMPTS[59]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 61) { ingate = bantest(si[61].tpratio,si[61].traceint,polyX513,polyY513,GNMPTS[60]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 62) { ingate = bantest(si[62].tpratio,si[62].traceint,polyX514,polyY514,GNMPTS[61]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 63) { ingate = bantest(si[63].tpratio,si[63].traceint,polyX515,polyY515,GNMPTS[62]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 64) { ingate = bantest(si[64].tpratio,si[64].traceint,polyX516,polyY516,GNMPTS[63]);
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

