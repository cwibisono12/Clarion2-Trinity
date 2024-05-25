//External Variables:
//C. Wibisono
//#include "global.h"


//////////////////////////////////////////
// INPUT CALIBRATION AND MAP PARAMETERS //
//////////////////////////////////////////

extern float ecal[][2];
extern float tcal[][2];

extern char map2type[];
extern int map2det[];
extern int map2deti[];
extern float map2beta[];
extern float mapangles[][2]; //theta and phi
extern float mapanglesi[][2]; //theta and phi
extern float mapangles1[][2]; //theta and phi
extern float mapangles2[][2];//theta and phi

//GAGG CALIBRATION:
extern int gaid[];
extern float gaggslope[];
extern float gaggintercept[];
extern float gaggquad[];
//INPUT BANANA GATE/////////////////////////////////////////////////////////////////////////////////////// //C. W
extern float polyX11[],polyY11[],polyX12[],polyY12[], 
polyX13[],polyY13[],polyX14[],polyY14[],
polyX15[],polyY15[],polyX16[],polyY16[],
polyX17[],polyY17[],polyX18[],polyY18[],
polyX21[],polyY21[],polyX22[],polyY22[],
polyX23[],polyY23[],polyX24[],polyY24[],
polyX25[],polyY25[],polyX26[],polyY26[],
polyX27[],polyY27[],polyX28[],polyY28[],
polyX29[],polyY29[],polyX210[],polyY210[],
polyX31[],polyY31[],polyX32[],polyY32[],
polyX33[],polyY33[],polyX34[],polyY34[],
polyX35[],polyY35[],polyX36[],polyY36[],
polyX37[],polyY37[],polyX38[],polyY38[],
polyX39[],polyY39[],polyX310[],polyY310[],
polyX311[],polyY311[],polyX312[],polyY312[],
polyX313[],polyY313[],polyX314[],polyY314[],
polyX41[],polyY41[],polyX42[],polyY42[],
polyX43[],polyY43[],polyX44[],polyY44[],
polyX45[],polyY45[],polyX46[],polyY46[],
polyX47[],polyY47[],polyX48[],polyY48[],
polyX49[],polyY49[],polyX410[],polyY410[],
polyX411[],polyY411[],polyX412[],polyY412[],
polyX413[],polyY413[],polyX414[],polyY414[],
polyX415[],polyY415[],polyX416[],polyY416[],
polyX51[],polyY51[],polyX52[],polyY52[],
polyX53[],polyY53[],polyX54[],polyY54[],
polyX55[],polyY55[],polyX56[],polyY56[],
polyX57[],polyY57[],polyX58[],polyY58[],
polyX59[],polyY59[],polyX510[],polyY510[],
polyX511[],polyY511[],polyX512[],polyY512[],
polyX513[],polyY513[],polyX514[],polyY514[],
polyX515[],polyY515[],polyX516[],polyY516[];
extern int GID[],GNMPTS[];
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//====================================================

//Kinematic Correction Parameters;
extern float Ebeam[];
extern float mbeam[];
extern float ml[];
extern float mh[];
extern float vbz[];

