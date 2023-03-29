#include <stdio.h>
#include "global.h"
//////////////////////////////////////////
// INPUT CALIBRATION AND MAP PARAMETERS //
//////////////////////////////////////////

float ecal[MAX_ID][2]={0};
float tcal[MAX_ID][2]={0};
float new_gain=1.0;

char map2type[MAX_ID]={0};
int map2det[MAX_ID]={0};
int map2deti[MAX_ID]={0};
float map2beta[MAX_ID]={0};
float mapangles[MAX_ID][2]={0}; //theta and phi
float mapanglesi[MAX_ID][2]={0}; //theta and phi
float mapangles1[MAX_ID][2]={0}; //theta and phi
float mapangles2[MAX_ID][2]={0};//theta and phi

//GAGG CALIBRATION:
int gaid[28]={0};
float gaggslope[28]={0};
float gaggintercept[28]={0};
float gaggquad[28]={0};
//INPUT BANANA GATE/////////////////////////////////////////////////////////////////////////////////////// //C. W
//INPUT BANANA GATE/////////////////////////////////////////////////////////////////////////////////////// //C. W
float polyX21[LINE_LENGTH]={0},polyY21[LINE_LENGTH]={0},polyX22[LINE_LENGTH]={0},polyY22[LINE_LENGTH]={0},
polyX23[LINE_LENGTH]={0},polyY23[LINE_LENGTH]={0},polyX24[LINE_LENGTH]={0},polyY24[LINE_LENGTH]={0},
polyX25[LINE_LENGTH]={0},polyY25[LINE_LENGTH]={0},polyX26[LINE_LENGTH]={0},polyY26[LINE_LENGTH]={0},
polyX27[LINE_LENGTH]={0},polyY27[LINE_LENGTH]={0},polyX28[LINE_LENGTH]={0},polyY28[LINE_LENGTH]={0},
polyX29[LINE_LENGTH]={0},polyY29[LINE_LENGTH]={0},polyX210[LINE_LENGTH]={0},polyY210[LINE_LENGTH]={0},
polyX41[LINE_LENGTH]={0},polyY41[LINE_LENGTH]={0},polyX42[LINE_LENGTH]={0},polyY42[LINE_LENGTH]={0},
polyX43[LINE_LENGTH]={0},polyY43[LINE_LENGTH]={0},polyX44[LINE_LENGTH]={0},polyY44[LINE_LENGTH]={0},
polyX45[LINE_LENGTH]={0},polyY45[LINE_LENGTH]={0},polyX46[LINE_LENGTH]={0},polyY46[LINE_LENGTH]={0},
polyX47[LINE_LENGTH]={0},polyY47[LINE_LENGTH]={0},polyX48[LINE_LENGTH]={0},polyY48[LINE_LENGTH]={0},
polyX49[LINE_LENGTH]={0},polyY49[LINE_LENGTH]={0},polyX410[LINE_LENGTH]={0},polyY410[LINE_LENGTH]={0},
polyX411[LINE_LENGTH]={0},polyY411[LINE_LENGTH]={0},polyX412[LINE_LENGTH]={0},polyY412[LINE_LENGTH]={0},
polyX413[LINE_LENGTH]={0},polyY413[LINE_LENGTH]={0},polyX414[LINE_LENGTH]={0},polyY414[LINE_LENGTH]={0},
polyX415[LINE_LENGTH]={0},polyY415[LINE_LENGTH]={0},polyX416[LINE_LENGTH]={0},polyY416[LINE_LENGTH]={0};
int GID[27]={0},GNMPTS[LINE_LENGTH]={0};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Charged Particle Type Parameter:
float Ebeam[5],mbeam[5],ml[5],mh[5],vbz[5];
