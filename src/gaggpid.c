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
    && (si[i].energy > 0. && si[i].energy < 4096.)
    && gemvalid > 0
	 ){
    if (i == 1){pid_qdc21[(int)si[1].energy][(int)si[i].tpratio]++;} 
    if (i == 2){pid_qdc22[(int)si[2].energy][(int)si[i].tpratio]++;} 
    if (i == 3){pid_qdc23[(int)si[3].energy][(int)si[i].tpratio]++;} 
    if (i == 4){pid_qdc24[(int)si[4].energy][(int)si[i].tpratio]++;} 
    if (i == 5){pid_qdc25[(int)si[5].energy][(int)si[i].tpratio]++;} 
    if (i == 6){pid_qdc26[(int)si[6].energy][(int)si[i].tpratio]++;} 
    if (i == 7){pid_qdc27[(int)si[7].energy][(int)si[i].tpratio]++;} 
    if (i == 8){pid_qdc28[(int)si[8].energy][(int)si[i].tpratio]++;} 
    if (i == 9){pid_qdc29[(int)si[9].energy][(int)si[i].tpratio]++;} 
    if (i == 10){pid_qdc210[(int)si[10].energy][(int)si[i].tpratio]++;} 
    if (i == 11){pid_qdc41[(int)si[11].energy][(int)si[i].tpratio]++;} 
    if (i == 12){pid_qdc42[(int)si[12].energy][(int)si[i].tpratio]++;} 
    if (i == 13){pid_qdc43[(int)si[13].energy][(int)si[i].tpratio]++;} 
    if (i == 14){pid_qdc44[(int)si[14].energy][(int)si[i].tpratio]++;} 
    if (i == 15){pid_qdc45[(int)si[15].energy][(int)si[i].tpratio]++;} 
    if (i == 16){pid_qdc46[(int)si[16].energy][(int)si[i].tpratio]++;} 
    if (i == 17){pid_qdc47[(int)si[17].energy][(int)si[i].tpratio]++;} 
    if (i == 18){pid_qdc48[(int)si[18].energy][(int)si[i].tpratio]++;} 
    if (i == 19){pid_qdc49[(int)si[19].energy][(int)si[i].tpratio]++;} 
    if (i == 20){pid_qdc410[(int)si[20].energy][(int)si[i].tpratio]++;} 
    if (i == 21){pid_qdc411[(int)si[21].energy][(int)si[i].tpratio]++;} 
    if (i == 22){pid_qdc412[(int)si[22].energy][(int)si[i].tpratio]++;} 
    if (i == 23){pid_qdc413[(int)si[23].energy][(int)si[i].tpratio]++;} 
    if (i == 24){pid_qdc414[(int)si[24].energy][(int)si[i].tpratio]++;} 
    if (i == 25){pid_qdc415[(int)si[25].energy][(int)si[i].tpratio]++;} 
    if (i == 26){pid_qdc416[(int)si[26].energy][(int)si[i].tpratio]++;} 
}

		}// end Si condition for GAGG firing:
    
	}//end Si

}

