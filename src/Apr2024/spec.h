///////////////////////////////////////
// SPECTRA and FILE NAME DEFINITIONS // 
///////////////////////////////////////
//All spectra are considered two-dimensional arrays
//Must add "write spectra" at end of file 

//[Y-dim][X-dim]

////////////////
//Event Spectra
////////////////

#define HIT "hit.spn"

#define MULT "mult.spn" //total detector multiplicity for one event

#define TDIFID "tdif.sec" //time diference between sequential events of a single channel ; 1 micro-second bins

#define E_RAW "e_raw.sec"

#define E_CAL "e_cal.sec"

#define TEVT_RAW "tevt_raw.sec" // 10 second bins

#define TEVT_CAL "tevt_cal.sec"

#define TCFD_RAW "tcfd_raw.sec"

#define TDIF_RAW "tdif_raw.spn" // 10 nano-second bins

#define TDIF_CAL "tdif_cal.spn" // 10 nano-second bins

#define TDIF_CAL0_ETHRESH "tdif_cal0_ethresh.spn" //time difference relative to channel 0 ; 10 nano-second bins

//#define IDID "idid.spn" //id vs id correlation hit pattern 


////////////////////////////
//Detector Processed Spectra
////////////////////////////

//Ge

#define GE_BGO_TDIF "ge_bgo_tdif.spn"

#define GE_XTL_TDIF "ge_xtl_tdif.spn"

#define GE_XTL_TDIF_ETHRESH "ge_xtl_tdif_ethresh.spn"

#define GE_SPE_XTL "ge_spe_xtl.sec"

#define GE_SPE "ge_spe.sec"

#define GE_SPE_DOPP "ge_spe_dopp.sec" //C. W

#define GE_SPE_CLEAN "ge_spe_clean.spn"

//////////////////////
//Final User Spectra
//////////////////////

//Singles-Gamma

#define PROTON_Gamma "p_gamm.spn2"

#define ALPHA_Gamma "a_gamm.spn2"

#define PROTON_ALPHA_Gamma "p_a_gamm.spn2"

#define Gamma "gamm.spn2"

//#define GAGG_GAMMA_TDIFF_RING2 "gagg_gamma_tdiff_ring2.spn"

//#define GAGG_GAMMA_TDIFF_RING4 "gagg_gamma_tdiff_ring4.spn"

#define GAGG_GAMMA_TDIFF2 "gagg_gamma_tdiff2.spn2"

#define GATED_GAMMA "gated_gamma.spn2"

#define GATED_GAMMATRUE "gated_gammatrue.spn2"

#define GATED_GAMMA_NP "gated_gamma_np.spn2"

#define GATED_GAMMA_NPTRUE "gated_gamma_nptrue.spn2"

#define E_LIT2 "e_lit2.spn2"

#define E_LIT4 "e_lit4.spn2"

#define THETA_RES2 "theta_res2.spn2"

#define THETA_RES4 "theta_res4.spn2"

#define GAGGHIST "gagghist.spn"
