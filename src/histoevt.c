#include <stdio.h>
#include "global.h"
#include "histoevt.h"


void pidevt(struct subevent *subevt){

float peak=(subevt->qsum[3] - (20./(31.+29.))*(subevt->qsum[0]+subevt->qsum[1]))/10.;
float tail=(subevt->qsum[5] - (55./(31.+29.))*(subevt->qsum[0]+subevt->qsum[1]))/10.;

float traceint=(subevt->qsum[2] + subevt->qsum[3] + subevt->qsum[4] + subevt->qsum[5] + subevt->qsum[6] - (115./(31.+29.))*(subevt->qsum[0] + subevt->qsum[1]))/30.;

float tpratio=4000.*tail/peak;

if (subevt->id == 60 || subevt->id ==70 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc21[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 61 || subevt->id ==71 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc22[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 62 || subevt->id ==72 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc23[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 63 || subevt->id ==73 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc24[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 64 || subevt->id ==74 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc25[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 65 || subevt->id ==75 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc26[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 66 || subevt->id ==76 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc27[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 67 || subevt->id ==77 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc28[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 68 || subevt->id ==78 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc29[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 69 || subevt->id ==79 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc210[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 80 || subevt->id == 96){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc41[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 81 || subevt->id == 97){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc42[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 82 || subevt->id == 98){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc43[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 83 || subevt->id == 99){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc44[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 84 || subevt->id == 100){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc45[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 85 || subevt->id == 101 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc46[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 86 || subevt->id == 102 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc47[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 87 || subevt->id == 103 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc48[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 88 || subevt->id == 104 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc49[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 89 || subevt->id == 105 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc410[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 90 || subevt->id == 106 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc411[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 91 || subevt->id == 107 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc412[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 92 || subevt->id == 108 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc413[(int) traceint][(int) tpratio]++;
}

if (subevt->id == 93 || subevt->id == 109 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc414[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 94 || subevt->id == 110 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc415[(int) traceint][(int) tpratio]++;
}


if (subevt->id == 95 || subevt->id == 111 ){
if((traceint > 0 && traceint < 4096) && (tpratio > 0 && tpratio <4096))
pid_qdc416[(int) traceint][(int) tpratio]++;
}



}
