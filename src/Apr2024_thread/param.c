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
int gaid[MAX_SI]={0};
float gaggslope[MAX_SI]={0};
float gaggintercept[MAX_SI]={0};
float gaggquad[MAX_SI]={0};
//INPUT BANANA GATE/////////////////////////////////////////////////////////////////////////////////////// //C. W
//INPUT BANANA GATE/////////////////////////////////////////////////////////////////////////////////////// //C. W
float polyX11[LINE_LENGTH]={0},polyY11[LINE_LENGTH]={0},polyX12[LINE_LENGTH]={0},polyY12[LINE_LENGTH]={0},
polyX13[LINE_LENGTH]={0},polyY13[LINE_LENGTH]={0},polyX14[LINE_LENGTH]={0},polyY14[LINE_LENGTH]={0},
polyX15[LINE_LENGTH]={0},polyY15[LINE_LENGTH]={0},polyX16[LINE_LENGTH]={0},polyY16[LINE_LENGTH]={0},
polyX17[LINE_LENGTH]={0},polyY17[LINE_LENGTH]={0},polyX18[LINE_LENGTH]={0},polyY18[LINE_LENGTH]={0},
polyX21[LINE_LENGTH]={0},polyY21[LINE_LENGTH]={0},polyX22[LINE_LENGTH]={0},polyY22[LINE_LENGTH]={0},
polyX23[LINE_LENGTH]={0},polyY23[LINE_LENGTH]={0},polyX24[LINE_LENGTH]={0},polyY24[LINE_LENGTH]={0},
polyX25[LINE_LENGTH]={0},polyY25[LINE_LENGTH]={0},polyX26[LINE_LENGTH]={0},polyY26[LINE_LENGTH]={0},
polyX27[LINE_LENGTH]={0},polyY27[LINE_LENGTH]={0},polyX28[LINE_LENGTH]={0},polyY28[LINE_LENGTH]={0},
polyX29[LINE_LENGTH]={0},polyY29[LINE_LENGTH]={0},polyX210[LINE_LENGTH]={0},polyY210[LINE_LENGTH]={0},
polyX31[LINE_LENGTH]={0},polyY31[LINE_LENGTH]={0},polyX32[LINE_LENGTH]={0},polyY32[LINE_LENGTH]={0},
polyX33[LINE_LENGTH]={0},polyY33[LINE_LENGTH]={0},polyX34[LINE_LENGTH]={0},polyY34[LINE_LENGTH]={0},
polyX35[LINE_LENGTH]={0},polyY35[LINE_LENGTH]={0},polyX36[LINE_LENGTH]={0},polyY36[LINE_LENGTH]={0},
polyX37[LINE_LENGTH]={0},polyY37[LINE_LENGTH]={0},polyX38[LINE_LENGTH]={0},polyY38[LINE_LENGTH]={0},
polyX39[LINE_LENGTH]={0},polyY39[LINE_LENGTH]={0},polyX310[LINE_LENGTH]={0},polyY310[LINE_LENGTH]={0},
polyX311[LINE_LENGTH]={0},polyY311[LINE_LENGTH]={0},polyX312[LINE_LENGTH]={0},polyY312[LINE_LENGTH]={0},
polyX313[LINE_LENGTH]={0},polyY313[LINE_LENGTH]={0},polyX314[LINE_LENGTH]={0},polyY314[LINE_LENGTH]={0},
polyX41[LINE_LENGTH]={0},polyY41[LINE_LENGTH]={0},polyX42[LINE_LENGTH]={0},polyY42[LINE_LENGTH]={0},
polyX43[LINE_LENGTH]={0},polyY43[LINE_LENGTH]={0},polyX44[LINE_LENGTH]={0},polyY44[LINE_LENGTH]={0},
polyX45[LINE_LENGTH]={0},polyY45[LINE_LENGTH]={0},polyX46[LINE_LENGTH]={0},polyY46[LINE_LENGTH]={0},
polyX47[LINE_LENGTH]={0},polyY47[LINE_LENGTH]={0},polyX48[LINE_LENGTH]={0},polyY48[LINE_LENGTH]={0},
polyX49[LINE_LENGTH]={0},polyY49[LINE_LENGTH]={0},polyX410[LINE_LENGTH]={0},polyY410[LINE_LENGTH]={0},
polyX411[LINE_LENGTH]={0},polyY411[LINE_LENGTH]={0},polyX412[LINE_LENGTH]={0},polyY412[LINE_LENGTH]={0},
polyX413[LINE_LENGTH]={0},polyY413[LINE_LENGTH]={0},polyX414[LINE_LENGTH]={0},polyY414[LINE_LENGTH]={0},
polyX415[LINE_LENGTH]={0},polyY415[LINE_LENGTH]={0},polyX416[LINE_LENGTH]={0},polyY416[LINE_LENGTH]={0},
polyX51[LINE_LENGTH]={0},polyY51[LINE_LENGTH]={0},polyX52[LINE_LENGTH]={0},polyY52[LINE_LENGTH]={0},
polyX53[LINE_LENGTH]={0},polyY53[LINE_LENGTH]={0},polyX54[LINE_LENGTH]={0},polyY54[LINE_LENGTH]={0},
polyX55[LINE_LENGTH]={0},polyY55[LINE_LENGTH]={0},polyX56[LINE_LENGTH]={0},polyY56[LINE_LENGTH]={0},
polyX57[LINE_LENGTH]={0},polyY57[LINE_LENGTH]={0},polyX58[LINE_LENGTH]={0},polyY58[LINE_LENGTH]={0},
polyX59[LINE_LENGTH]={0},polyY59[LINE_LENGTH]={0},polyX510[LINE_LENGTH]={0},polyY510[LINE_LENGTH]={0},
polyX511[LINE_LENGTH]={0},polyY511[LINE_LENGTH]={0},polyX512[LINE_LENGTH]={0},polyY512[LINE_LENGTH]={0},
polyX513[LINE_LENGTH]={0},polyY513[LINE_LENGTH]={0},polyX514[LINE_LENGTH]={0},polyY514[LINE_LENGTH]={0},
polyX515[LINE_LENGTH]={0},polyY515[LINE_LENGTH]={0},polyX516[LINE_LENGTH]={0},polyY516[LINE_LENGTH]={0};
int GID[MAX_SI]={0},GNMPTS[LINE_LENGTH]={0};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Charged Particle Type Parameter:
float Ebeam[5],mbeam[5],ml[5],mh[5],vbz[5];
