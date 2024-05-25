#include <stdio.h>

extern unsigned long long int pileupcount;

extern float new_gain;

extern float ecal[][2];

extern float tcal[][2];

///////////////////////////////////////
// SPECTRA and FILE NAME DEFINITIONS // 
///////////////////////////////////////
//All spectra are considered two-dimensional arrays
//Must add "write spectra" at end of file 

//[Y-dim][X-dim]

////////////////
//Event Spectra
////////////////

extern int hit[][4096]; //first for all hits, second for pilup hits

extern int mult[][4096];

extern int tdifid[][8192];

extern int e_raw[][8192];

extern int e_cal[][8192];

extern int tevt_raw[][8192];

extern int tevt_cal[][8192]; // 10 second bins

extern int tcfd_raw[][8192];

extern int tdif_raw[][4096];

extern int tdif_cal[][4096];

extern int tdif_cal0_ethresh[][4096];

//extern int idid[][];


////////////////////////////
//Detector Processed Spectra
////////////////////////////

//Ge

extern int ge_bgo_tdif[][4096];

extern int ge_xtl_tdif[][4096];

extern int ge_xtl_tdif_ethresh[][4096];

extern int ge_spe_xtl[][8192];

extern int ge_spe[][8192];

extern int ge_spe_dopp[][8192];

extern int ge_spe_clean[][8192];

//int pid_all[26][4096][4096]={{{0}}}; //TJG + CW
//////////////////////
//Final User Spectra
//////////////////////

//Singles-Gamma

extern int p_gamm[][8192];

extern int a_gamm[][8192];

extern int p_a_gamm[][8192];

extern int gamm[][8192];

//GAGG-GAMMA
//extern int gagg_gamma_tdiff_ring2[][4096];
//extern int gagg_gamma_tdiff_ring4[][4096];
extern int gagg_gamma_tdiff2[][5000];
extern int gated_gamma[][5000];
extern int gated_gammatrue[][5000];
extern int gated_gamma_np[][5000];
extern int gated_gamma_nptrue[][5000];
extern int e_lit2[][4096];
extern int e_lit4[][4096];
extern int theta_res2[][4096];
extern int theta_res4[][4096];
extern int gagghist[][4096];

extern int gated_gammaRing[][5000];
extern int gated_gammaRingtrue[][5000];
extern int gated_gammaRing_np[][5000];
extern int gated_gammaRing_nptrue[][5000];

