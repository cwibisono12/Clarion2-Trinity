#include <stdio.h>
#include "global.h"
///////////////////////////////////////
// SPECTRA and FILE NAME DEFINITIONS // 
///////////////////////////////////////
//All spectra are considered two-dimensional arrays
//Must add "write spectra" at end of file 

//[Y-dim][X-dim]

////////////////
//Event Spectra
////////////////

int hit[2][4096]={0}; //first for all hits, second for pilup hits

int mult[1][4096]={0};

int tdifid[MAX_ID][8192]={0};

int e_raw[MAX_ID][8192]={0};

int e_cal[MAX_ID][8192]={0};

int tevt_raw[MAX_ID][8192]={0};

int tevt_cal[MAX_ID][8192]={0}; // 10 second bins

int tcfd_raw[MAX_ID][8192]={0};

int tdif_raw[MAX_ID][4096]={0};

int tdif_cal[MAX_ID][4096]={0};

int tdif_cal0_ethresh[MAX_ID][4096]={0};

//int idid[MAX_ID][MAX_ID]={0};


////////////////////////////
//Detector Processed Spectra
////////////////////////////

//Ge

int ge_bgo_tdif[MAX_GE][4096]={0};

int ge_xtl_tdif[MAX_GE][4096]={0};

int ge_xtl_tdif_ethresh[MAX_GE][4096]={0};

int ge_spe_xtl[MAX_GE*MAX_GE_XTL][8192]={0};

int ge_spe[MAX_GE][8192]={0};

int ge_spe_dopp[MAX_GE][8192]={0};

int ge_spe_clean[MAX_GE][8192]={0};

//int pid_all[26][4096][4096]={{{0}}}; //TJG + CW
//////////////////////
//Final User Spectra
//////////////////////

//Singles-Gamma

int p_gamm[10][8192]={0};

int a_gamm[10][8192]={0};

int p_a_gamm[1][8192]={0};

int gamm[1][8192]={0};


//GAGG-GAMMA
//int gagg_gamma_tdiff_ring2[4096][4096];
//int gagg_gamma_tdiff_ring4[4096][4096];
int gagg_gamma_tdiff2[5000][5000];
int gated_gamma[16][5000];
int gated_gammatrue[16][5000];
int gated_gamma_np[16][5000];
int gated_gamma_nptrue[16][5000];
int e_lit2[MAX_SI][4096];
int e_lit4[MAX_SI][4096];
int theta_res2[MAX_SI][4096];
int theta_res4[MAX_SI][4096];
int gagghist[1][4096];

int gated_gammaRing[MAX_SI][5000];
int gated_gammaRingtrue[MAX_SI][5000];
int gated_gammaRing_np[MAX_SI][5000];
int gated_gammaRing_nptrue[MAX_SI][5000];
