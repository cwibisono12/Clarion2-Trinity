//Usage: GAGG Data Processing: --Generate PID
//C. Wibisono
#include "global.h"
#include "mapping.h"
#include <math.h>
#include <stdbool.h>
#include "gaggpid.h"

void gaggpid(struct gdetector *ge, struct sidetector *si, int parttype,int gmult,int Egammin, int Egammax){
int gaggvalid=0;
int sicount;
int i,k;
int gemvalid=0;

//Checking whether an event is within a gamma gate window:
for(k=0;k<gmult;k++){
if(ge[k].energy >= Egammin && ge[k].energy <= Egammax)
gemvalid++;
}
//======================================================

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
if (parttype ==1) si[i].energy=100.*(gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
if (parttype ==2) si[i].energy=100.*(gaggquad[i-1]*pow(si[i].traceint,2.0)+gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
}


//SiPM angle assignment:
for (sicount=0;sicount<3;sicount++){
si[i].theta[sicount]=si[i].sitheta[0][sicount];
si[i].phi[sicount]=si[i].siphi[0][sicount];
//if (sicount == 0)
//printf("theta: %f phi: %f\n",si[i].theta[sicount],si[i].phi[sicount]);
}



//Generating PID:
if ((si[i].peak > 0. && si[i].peak < 4096.) && (si[i].tail > 0. && si[i].tail < 4096.) &&
    (si[i].traceint > 0. && si[i].traceint < 4096.) && (si[i].tpratio > 0. && si[i].tpratio < 4096.)
   
	 ){
    if (i == 1){pid_qdc11[(int)si[1].traceint][(int)si[i].tpratio]++;} 
    if (i == 2){pid_qdc12[(int)si[2].traceint][(int)si[i].tpratio]++;} 
    if (i == 3){pid_qdc13[(int)si[3].traceint][(int)si[i].tpratio]++;} 
    if (i == 4){pid_qdc14[(int)si[4].traceint][(int)si[i].tpratio]++;} 
    if (i == 5){pid_qdc15[(int)si[5].traceint][(int)si[i].tpratio]++;} 
    if (i == 6){pid_qdc16[(int)si[6].traceint][(int)si[i].tpratio]++;} 
    if (i == 7){pid_qdc17[(int)si[7].traceint][(int)si[i].tpratio]++;} 
    if (i == 8){pid_qdc18[(int)si[8].traceint][(int)si[i].tpratio]++;} 
    if (i == 9){pid_qdc21[(int)si[9].traceint][(int)si[i].tpratio]++;} 
    if (i == 10){pid_qdc22[(int)si[10].traceint][(int)si[i].tpratio]++;} 
    if (i == 11){pid_qdc23[(int)si[11].traceint][(int)si[i].tpratio]++;} 
    if (i == 12){pid_qdc24[(int)si[12].traceint][(int)si[i].tpratio]++;} 
    if (i == 13){pid_qdc25[(int)si[13].traceint][(int)si[i].tpratio]++;} 
    if (i == 14){pid_qdc26[(int)si[14].traceint][(int)si[i].tpratio]++;} 
    if (i == 15){pid_qdc27[(int)si[15].traceint][(int)si[i].tpratio]++;} 
    if (i == 16){pid_qdc28[(int)si[16].traceint][(int)si[i].tpratio]++;} 
    if (i == 17){pid_qdc29[(int)si[17].traceint][(int)si[i].tpratio]++;} 
    if (i == 18){pid_qdc210[(int)si[18].traceint][(int)si[i].tpratio]++;} 
    if (i == 19){pid_qdc31[(int)si[19].traceint][(int)si[i].tpratio]++;} 
    if (i == 20){pid_qdc32[(int)si[20].traceint][(int)si[i].tpratio]++;} 
    if (i == 21){pid_qdc33[(int)si[21].traceint][(int)si[i].tpratio]++;} 
    if (i == 22){pid_qdc34[(int)si[22].traceint][(int)si[i].tpratio]++;} 
    if (i == 23){pid_qdc35[(int)si[23].traceint][(int)si[i].tpratio]++;} 
    if (i == 24){pid_qdc36[(int)si[24].traceint][(int)si[i].tpratio]++;} 
    if (i == 25){pid_qdc37[(int)si[25].traceint][(int)si[i].tpratio]++;} 
    if (i == 26){pid_qdc38[(int)si[26].traceint][(int)si[i].tpratio]++;} 
    if (i == 27){pid_qdc39[(int)si[27].traceint][(int)si[i].tpratio]++;} 
    if (i == 28){pid_qdc310[(int)si[28].traceint][(int)si[i].tpratio]++;} 
    if (i == 29){pid_qdc311[(int)si[29].traceint][(int)si[i].tpratio]++;} 
    if (i == 30){pid_qdc312[(int)si[30].traceint][(int)si[i].tpratio]++;} 
    if (i == 31){pid_qdc313[(int)si[31].traceint][(int)si[i].tpratio]++;} 
    if (i == 32){pid_qdc314[(int)si[32].traceint][(int)si[i].tpratio]++;} 
    if (i == 33){pid_qdc41[(int)si[33].traceint][(int)si[i].tpratio]++;} 
    if (i == 34){pid_qdc42[(int)si[34].traceint][(int)si[i].tpratio]++;} 
    if (i == 35){pid_qdc43[(int)si[35].traceint][(int)si[i].tpratio]++;} 
    if (i == 36){pid_qdc44[(int)si[36].traceint][(int)si[i].tpratio]++;} 
    if (i == 37){pid_qdc45[(int)si[37].traceint][(int)si[i].tpratio]++;} 
    if (i == 38){pid_qdc46[(int)si[38].traceint][(int)si[i].tpratio]++;} 
    if (i == 39){pid_qdc47[(int)si[39].traceint][(int)si[i].tpratio]++;} 
    if (i == 40){pid_qdc48[(int)si[40].traceint][(int)si[i].tpratio]++;} 
    if (i == 41){pid_qdc49[(int)si[41].traceint][(int)si[i].tpratio]++;} 
    if (i == 42){pid_qdc410[(int)si[42].traceint][(int)si[i].tpratio]++;} 
    if (i == 43){pid_qdc411[(int)si[43].traceint][(int)si[i].tpratio]++;} 
    if (i == 44){pid_qdc412[(int)si[44].traceint][(int)si[i].tpratio]++;} 
    if (i == 45){pid_qdc413[(int)si[45].traceint][(int)si[i].tpratio]++;} 
    if (i == 46){pid_qdc414[(int)si[46].traceint][(int)si[i].tpratio]++;} 
    if (i == 47){pid_qdc415[(int)si[47].traceint][(int)si[i].tpratio]++;} 
    if (i == 48){pid_qdc416[(int)si[48].traceint][(int)si[i].tpratio]++;} 
    if (i == 49){pid_qdc51[(int)si[49].traceint][(int)si[i].tpratio]++;} 
    if (i == 50){pid_qdc52[(int)si[50].traceint][(int)si[i].tpratio]++;} 
    if (i == 51){pid_qdc53[(int)si[51].traceint][(int)si[i].tpratio]++;} 
    if (i == 52){pid_qdc54[(int)si[52].traceint][(int)si[i].tpratio]++;} 
    if (i == 53){pid_qdc55[(int)si[53].traceint][(int)si[i].tpratio]++;} 
    if (i == 54){pid_qdc56[(int)si[54].traceint][(int)si[i].tpratio]++;} 
    if (i == 55){pid_qdc57[(int)si[55].traceint][(int)si[i].tpratio]++;} 
    if (i == 56){pid_qdc58[(int)si[56].traceint][(int)si[i].tpratio]++;} 
    if (i == 57){pid_qdc59[(int)si[57].traceint][(int)si[i].tpratio]++;} 
    if (i == 58){pid_qdc510[(int)si[58].traceint][(int)si[i].tpratio]++;} 
    if (i == 59){pid_qdc511[(int)si[59].traceint][(int)si[i].tpratio]++;} 
    if (i == 60){pid_qdc512[(int)si[60].traceint][(int)si[i].tpratio]++;} 
    if (i == 61){pid_qdc513[(int)si[61].traceint][(int)si[i].tpratio]++;} 
    if (i == 62){pid_qdc514[(int)si[62].traceint][(int)si[i].tpratio]++;} 
    if (i == 63){pid_qdc515[(int)si[63].traceint][(int)si[i].tpratio]++;} 
    if (i == 64){pid_qdc516[(int)si[64].traceint][(int)si[i].tpratio]++;} 
}

		}// end Si condition for GAGG firing:
    
	}//end Si

}

