/**********************************************************/
/* PXI SCAN CODE -- J.M. Allmond (ORNL) -- July 2016      */
/* with some modifications-- C. Wibisono                  */
/* modifications include:                                 */
/* Doppler Corrections                                    */
/* Dopp-Nodopp &Dopp-Dopp Matrices                        */
/* PID for QDC sum                                        */
/* Incorporation of PID & Gamma                           */
/* !unpak data from Pixie-16 digitizers, event build,     */
/* !and create detctors and user defined spectra          */
/*                                                        */
/* gcc -o pxi-scan pxi-scan.c                             */
/* ./pxi-scan -op datafile calibrationfile mapfile        */
/*                                                        */
/* ..... calibration file optional                        */
/* ..... map file optional                                */
/* ..... u for update spectra                             */
/* ..... o for overwrite spectra                          */
/* ..... p for print realtime stats                       */
/**********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define PRINT_CAL 1
#define PRINT_MAP 1

#define DB(x) fwrite(x, sizeof(x), 1, debugfile);
#define DEBUGFN "debug.mat"

#define RAND ((float) rand() / ((unsigned int) RAND_MAX + 1))   // random number in interval (0,1)
#define TRUE  1
#define FALSE 0

#define LINE_LENGTH 120

#define MAX_CRATES 2
#define MAX_BOARDS_PER_CRATE 13
#define MAX_CHANNELS_PER_BOARD 16
#define BOARD_START 2

#define MAX_ID MAX_CRATES*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD

#define PIR 3.14159265/180.0 //C.W 

#define HEADER_LENGTH 4     //unit = words with 4 bytes per word
#define MAX_SUB_LENGTH 2016 //unit = words with 4 bytes per word ; 2004 --> 40 micro-second trace + 4 word header

#define EVENT_BUILD_TIME 190 // 100 = 1 micro-second ; should be < L + G ~ 5.04 us (note 0.08 us scale factor in set file)

#define RAWE_REBIN_FACTOR 2.0 // Rebin 32k pixie16 spectra to something smaller to fit better into 8k.

bool bantest(float x, float y, float *polyX, float *polyY, int polyCorners); //C.W
void write_data2(char *filename, short *data, int xdim, int ydim, int overwrite); //C.W
void write_data4(char *filename, int *data, int xdim, int ydim, int overwrite); //C.W

/////////////////////
// RAW EVENT TYPES //
/////////////////////
struct subevent
{
    int chn; 
    int sln;
    int crn;
    int id;
    int hlen;
    int elen;
    int trlen;          //number of samples
    int trwlen;         //number of words (two samples per word)
    int fcode;          //pileup flag
    long long int time;
    int ctime;
    int ctimef;
    int energy;
    int extra;
    short tr[4096];
    int esum[4];
    int qsum[8];        
}; 
struct subevent subevt[MAX_ID]={0};
int sevtmult=0;
unsigned long long int sevtcount=0;
unsigned long long int pileupcount=0;
unsigned long long int evtcount=0;

int plt200=0;
int pgt200=0;

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



////////////////////
// DETECTOR TYPES //
////////////////////

//G = Ge
#define MAX_GE 16 // max number of Ge detectors
#define MAX_GE_XTL 4 // max number of crystals per Ge detector
#define MAX_GE_SEG 3 // max number of segments per Ge detector
#define MAX_GE_BGO 1 // max number of BGO PMTs per Ge detector 
#define GE_BGO_SUPPRESSION TRUE
struct gdetector
{    
    int xmult;
    int xid[MAX_GE_XTL]; 
    int xe[MAX_GE_XTL];
    float xedopp[MAX_GE_XTL];   //C. W
    long long int xt[MAX_GE_XTL]; 
    int xct[MAX_GE_XTL];   
    float xtheta[MAX_GE_XTL][4]; 
    float xphi[MAX_GE_XTL][4];
    bool xpileup[MAX_GE_XTL];
    bool xsuppress[MAX_GE_XTL];
    int xsubevtid[MAX_GE_XTL];

    int smult;
    int sid[MAX_GE_SEG]; 
    int se[MAX_GE_SEG];
    long long int st[MAX_GE_SEG];
    int sct[MAX_GE_SEG];
    float stheta[MAX_GE_SEG][4]; 
    float sphi[MAX_GE_SEG][4];
    bool spileup[MAX_GE_SEG];
    bool ssuppress[MAX_GE_SEG];
    int ssubevtid[MAX_GE_SEG];

    int bgomult;
    int bgoid[MAX_GE_BGO];
    int bgoe[MAX_GE_BGO];
    long long int bgot[MAX_GE_BGO];
    int bgoct[MAX_GE_BGO];
    float bgotheta[MAX_GE_XTL][4]; 
    float bgophi[MAX_GE_XTL][4];
    bool bgopileup[MAX_GE_BGO];
    int bgosubevtid[MAX_GE_BGO];

    int id;
    int energy;
    int edop;
    long long int time;   
    int ctime;  
    float theta[3]; //det, xtl, or seg angle
    float phi[3]; //det, xtl, or seg angle  
    bool suppress; //at least one xtl was suppressed by bgo
    bool pileup; //two or more unspressed xtls but at least one had pileup
    bool nonprompt; //two or more unspressed xtls but at least one was non-prompt with first xtl
    bool clean;
    int validp;  //Flag for prompt //C.W
    int validnp; //Flag for nprompt //C.W
    float Dg[3]; //The unit vector of the target to Ge Detector
    float etrue; //true energy after doppler correct
    float etrue2;
    float betavalue;
    float alpha; //cos(part-gamma)

//attribute for polarization analysis:
    int xvalid;
    int xidvalid[2];
    int xevalid[2];
   // float thetavalid[2];
    int epol;
    int epoldopp;
}; 
struct gdetector ge[MAX_GE]={0};
int gmult=0;
unsigned long long int gcount=0;

//S = Si============================       //C. W
#define MAX_SI_PAIR 2
#define MAX_SI 28
struct sidetector
{    
    int simult;
    int siid[MAX_SI_PAIR]; 
    int sie[MAX_SI_PAIR];
    long long int sit[MAX_SI_PAIR]; 
    int sict[MAX_SI_PAIR];   
    float sitheta[MAX_SI_PAIR][4]; 
    float siphi[MAX_SI_PAIR][4];
    bool sipileup[MAX_SI_PAIR];
    bool sisuppress[MAX_SI_PAIR];
    int sisubevtid[MAX_SI_PAIR];
    int siqdc[MAX_SI_PAIR][8];   

    int id;
    //int energy;
    float energy; //initially integer
    float peak;
    float tail;
    long long int time;   
    int ctime;  
    float theta[3]; //det, xtl, or seg angle
    float phi[3]; //det, xtl, or seg angle  
    bool suppress; //at least one xtl was suppressed by bgo
    bool pileup; //two or more unspressed xtls but at least one had pileup
    bool nonprompt; //two or more unspressed xtls but at least one was non-prompt with first xtl
    bool clean;
    int valid;
   // int traceint;
   float traceint; //initially integer 
   float tpratio; //tail to peak ratio
   float vl[3]; //velocity of light particle
   float vh[3]; //velocity of heavy particle
   float betarest; //betavalue of residual particle
   float thetarest; //thetaresidual in degree
}; 
struct sidetector si[MAX_SI]={0};
//int simult=0;
unsigned long long int sicount=0;
//====================================================



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
int hit[2][4096]={0}; //first for all hits, second for pilup hits

#define MULT "mult.spn" //total detector multiplicity for one event
int mult[1][4096]={0};

#define TDIFID "tdif.sec" //time diference between sequential events of a single channel ; 1 micro-second bins
int tdifid[MAX_ID][8192]={0};

#define E_RAW "e_raw.sec"
int e_raw[MAX_ID][8192]={0};

#define E_CAL "e_cal.sec"
int e_cal[MAX_ID][8192]={0};

#define TEVT_RAW "tevt_raw.sec" // 10 second bins
int tevt_raw[MAX_ID][8192]={0};

#define TEVT_CAL "tevt_cal.sec"
int tevt_cal[MAX_ID][8192]={0}; // 10 second bins

#define TCFD_RAW "tcfd_raw.sec"
int tcfd_raw[MAX_ID][8192]={0};

#define TDIF_RAW "tdif_raw.spn" // 10 nano-second bins
int tdif_raw[MAX_ID][4096]={0};

#define TDIF_CAL "tdif_cal.spn" // 10 nano-second bins
int tdif_cal[MAX_ID][4096]={0};

#define TDIF_CAL0_ETHRESH "tdif_cal0_ethresh.spn" //time difference relative to channel 0 ; 10 nano-second bins
int tdif_cal0_ethresh[MAX_ID][4096]={0};

#define IDID "idid.spn" //id vs id correlation hit pattern 
int idid[MAX_ID][MAX_ID]={0};


////////////////////////////
//Detector Processed Spectra
////////////////////////////

//Ge

#define GE_BGO_TDIF "ge_bgo_tdif.spn"
int ge_bgo_tdif[MAX_GE][4096]={0};

#define GE_XTL_TDIF "ge_xtl_tdif.spn"
int ge_xtl_tdif[MAX_GE][4096]={0};

#define GE_XTL_TDIF_ETHRESH "ge_xtl_tdif_ethresh.spn"
int ge_xtl_tdif_ethresh[MAX_GE][4096]={0};

#define GE_SPE_XTL "ge_spe_xtl.sec"
int ge_spe_xtl[MAX_GE*MAX_GE_XTL][8192]={0};

#define GE_SPE "ge_spe.sec"
int ge_spe[MAX_GE][8192]={0};

#define GE_SPE_DOPP "ge_spe_dopp.sec" //C. W
int ge_spe_dopp[MAX_GE][8192]={0};

#define GE_SPE_CLEAN "ge_spe_clean.spn"
int ge_spe_clean[MAX_GE][4096]={0};

//int pid_all[26][4096][4096]={{{0}}}; //TJG + CW
//////////////////////
//Final User Spectra
//////////////////////

//gamma-gamma

/*

#define GG_TDIF "gg_tdif.spn"
int gg_tdif[4096][4096]={0};

#define GG_TDIF1 "gg_tdif1.spn"
int gg_tdif1[4096][4096]={0};

#define GG_TDIF2 "gg_tdif2.spn"
int gg_tdif2[4096][4096]={0};

#define GG_PROMPT "gg_prompt.spn2"
int gg_prompt[5000][5000]={0};

#define GG_PROMPT_PURE "gg_prompt_pure.spn2"
int gg_prompt_pure[5000][5000]={0}; 

#define GG_NPROMPT_PURE "gg_nprompt_pure.spn2"
int gg_nprompt_pure[5000][5000]={0};

#define GG_PROMPT_nddopp "gg_prompt_nddopp.spn2" 
int gg_prompt_nddopp[5000][5000]={0};  //C.W

#define GG_PROMPT_ddopp "gg_prompt_ddopp.spn2"
int gg_prompt_ddopp[5000][5000]={0}; //C.W

#define GG_NPROMPT "gg_nprompt.spn2"
int gg_nprompt[5000][5000]={0};


#define GG_NPROMPT1 "gg_nprompt1.spn2"
int gg_nprompt1[5000][5000]={0};

#define GG_NPROMPT2 "gg_nprompt2.spn2"
int gg_nprompt2[5000][5000]={0};

#define GG_NPROMPT_nddopp "gg_nprompt_nddopp.spn2"
int gg_nprompt_nddopp[5000][5000]={0}; //C.W

#define GG_NPROMPT_nddopp1 "gg_nprompt_nddopp1.spn2"
int gg_nprompt_nddopp1[5000][5000]={0}; //C.W

#define GG_NPROMPT_ddopp1 "gg_nprompt_ddop1.spn2"
int gg_nprompt_ddopp1[5000][5000]={0}; //C.W

#define GG_NPROMPT_nddopp2 "gg_nprompt_nddopp2.spn2"
int gg_nprompt_nddopp2[5000][5000]={0}; //C.W

#define GG_NPROMPT_ddopp2 "gg_nprompt_ddopp2.spn2"
int gg_nprompt_ddopp2[5000][5000]={0}; //C.W

#define GG_NPROMPT_ddopp "gg_nprompt_ddopp.spn2"
int gg_nprompt_ddopp[5000][5000]={0}; //C.W

*/
#define PID_Gammaa "gagg_gamma_tdiff_ring2.spn"
int gagg_gamma_tdiff_ring2[4096][4096]={0}; //C.W + T.J 

#define PID_Gammab "gagg_gamma_tdiff_ring4.spn"
int gagg_gamma_tdiff_ring4[4096][4096]={0};

#define PID_Gamma2 "gagg_gamma_tdiff2.spn2"
int gagg_gamma_tdiff2[5000][5000]={0};

#define PID_Gamma_R "ge_gagg_spe_r.spn"
int ge_gagg_spe_r[MAX_GE][4096]={0}; //C.W

#define PID_Gamma_1Dp "gated_gamma.spn2"
int gated_gamma[1][5000]={0}; //C. W + T.J

#define PID_Gamma_1Dptrue "gated_gammatrue.spn2"
int gated_gammatrue[1][5000]={0};

#define PID_Gamma_1Dnp "gated_gamma_np.spn2"
int gated_gamma_np[1][5000]={0}; //C.W + T.J

#define PID_Gamma_1Dnptrue "gated_gamma_nptrue.spn2"
int gated_gamma_nptrue[1][5000]={0};


#define GAGG_HIST "gagghist.spn"
int gagghist[1][4096]={0};

#define THETA_RES2 "theta_res2.spn"
int theta_res2[MAX_SI][4096]={0};

#define E_LIT2 "e_lit2.spn"
int e_lit2[MAX_SI][4096]={0};

#define E_LIT4 "e_lit4.spn"
int e_lit4[MAX_SI][4096]={0};

#define THETA_RES4 "theta_res4.spn"
int theta_res4[MAX_SI][4096]={0};

#define POL_PARA "pol_para.spn2"
int pol_para[5000][5000]={0};

#define POL_PERP "pol_perp.spn2"
int pol_perp[5000][5000]={0};



void detmaps(int sevtmult, struct subevent *subevt, struct gdetector *ge, struct sidetector *si);
int gaggproc(struct sidetector *si);
int gamproc(struct gdetector *ge, int ge_spe_xtl[MAX_GE_XTL*MAX_GE][8192], int ge_spe[MAX_GE][8192]);
void doppcorrect(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int option, float ml, float mbeam, float mh, float vbz, float betascale);
void pol(int gaggvalid, int gmult, struct gdetector *ge, int pol_para[5000][5000], int pol_perp[5000][5000]);
//void doppcorrect(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int option, int betascale);
void gammagagg1Dspectra(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int gagg_gamma_tdiff_ring2[4096][4096], int gagg_gamma_tdiff_ring4[4096][4096], int gagg_gamma_tdiff2[5000][5000],int gated_gamma[1][5000], int gated_gammatrue[1][5000], int gated_gamma_np[1][5000], int gated_gamma_nptrue[1][5000]);
void parthit(int gaggvalid, int gmult, struct gdetector *ge, struct sidetector *si, int e_lit2[MAX_SI][4096], int e_lit4[MAX_SI][4096], int theta_res2[MAX_SI][4096], int theta_res4[MAX_SI][4096], int gagghist[1][4096]);

///////////////////////////////////
// START OF MAIN FUNCTION        //
///////////////////////////////////
int main(int argc, char **argv) {
  
  int i=0, j=0, k=0, l=0; //l is for iteration over trace
  float tempf=0;
  int max1=0, min1=0;
  int max2=0, min2=0;
  int maxid1=-1, minid1=-1;
  int maxid2=-1, minid2=-1;
  div_t e_div;
  lldiv_t lle_div;
//Parameter of kinematics Correction:
  float Ebeam,mbeam,ml,mh,amu;
  amu=931.5;
  Ebeam=30;
  mbeam=16*amu;
  ml=1*amu;
  mh=33*amu;
  float vbz=pow(2*Ebeam/mbeam,0.5);
  float betascale=atof(argv[7]);
  int option=atoi(argv[8]); //use 1 to use doppler corection with kinmat or 2 to use old doppler
//  int Emin=atoi(argv[6]); //C.W
//  int Emax=atof(argv[7]); //C.W
  int overwrite = 1;    
 

  double etrace0,etrace1,btrace0,btrace1;
  double ptrace0,ptrace1,ttrace0,ttrace1,tautrace0,tautrace1;
  int dbcount = 0;
  long long int strace[500];
  memset(strace, 0, sizeof(strace));

  //temp buffer for each sub event
  unsigned int sub[MAX_SUB_LENGTH];
  memset(sub, 0, sizeof(sub));
  
  //Reference time and difference for event building
  long long int etime, tdif, idtime[MAX_ID]={0}, temptime;
  
   
  // Check that the corrent number of arguments were provided.
  if (argc<3)    {
    printf("Incorrect number of arguments:\n");
    printf("%s -op datafile calibrationfile mapfile \n", argv[0]);
    printf("\n .... calibration file is optional\n");
    printf(" .... map file is optional\n");
    printf(" .... o for overwrite spectra\n");
    printf(" .... u for update spectra\n");
    printf(" .... p for print realtime stats\n");
    return 1;
  }
    
  if(strstr(argv[1], "u") != NULL) {
    overwrite = 0;
    printf("Updating Spectra\n");
  }
  else {
    printf("Overwriting Spectra\n");
  }  
  
  //open list-mode data file from PXI digitizer  
  FILE *fpr;
  long int fprsize,fprpos;
  if ((fpr = fopen(argv[2], "r")) == NULL) {
    fprintf(stderr, "Error, cannot open input file %s\n", argv[2]);
    return 1;
  }
 
  //get file size
  fseek(fpr, 0L, SEEK_END);
  fprsize = ftell(fpr);
  rewind(fpr);
  
   
  //open debug file for streaming an 1d array
  FILE *debugfile;
  if ((debugfile = fopen(DEBUGFN, "w")) == NULL) {
    fprintf(stderr, "Error, cannot open %s\n", DEBUGFN);
    return 1;
  }   
   

  //buffer for reading in text files for calibrations and maps below
  char line[LINE_LENGTH];


  //open energy and time calibration file (e.g., *.ca3 file)  
  int calid=0; 
  float caloffset=0.0, calgain=0.0;
  int firstcal=0, necal=0, ntcal=0;
  
  FILE *fprcal;
  
  for (i=0; i<MAX_ID; i++) {
    ecal[i][0] = caloffset;
    ecal[i][1] = calgain;
    tcal[i][0] = caloffset;
    tcal[i][1] = calgain;     
  }  
  
  if (argc >= 4) {

    if ((fprcal = fopen(argv[3], "r")) == NULL) {
        fprintf(stderr, "Error, cannot open input file %s\n", argv[3]);
        return 1;
    } 
    
    printf("%s loaded!\n", argv[3]);

	while(fgets(line, LINE_LENGTH, fprcal) != NULL){
	    calid=0; caloffset=0; calgain=0;
        for(i=0; i<LINE_LENGTH; i++){
            if(line[i]=='#'){
                if(PRINT_CAL)printf("%s", line);
                break;
            }
            else if(line[i]>=0){
                if (firstcal==0) {
                    sscanf(line,"%f\n", &new_gain);    
	                if(PRINT_CAL) printf("%f\n", new_gain);                                       
                    firstcal=1;
                    break;
                }    
	            sscanf(line,"%d\t%f\t%f\n", &calid, &caloffset, &calgain);
	            if(PRINT_CAL) printf("%d\t%.4f\t%.4f\n", calid, caloffset, calgain);
		        if(calid >=0 && calid < MAX_ID) {
                    ecal[calid][0] = caloffset;
                    ecal[calid][1] = calgain;
                    necal++;
          		    break;
			    }
		        if(calid >=1000 && calid < 1000+MAX_ID) {
                    tcal[calid-1000][0] = caloffset;
                    tcal[calid-1000][1] = calgain;
                    ntcal++;
          		    break;
			    }			    
			    else {
				    printf("Error in reading %s : bad id or format\n", argv[3]);
				    return -1;
			    }
            }
	        else if(line[i]=='\n'){
                if(PRINT_CAL) printf("\n");
                break;
            }
            else {
                continue;
            }

        }
        	memset(line, 0, LINE_LENGTH);
 	  }
 	  
 	fclose(fprcal);
    printf("read %d energy calibrations\n", necal);
    printf("read %d time calibrations\n", ntcal);
  }
  


  //open ID->Detector map file (e.g., *.map file)  
  int mapid=0, detid=0, detidi; 
  char dettype=0;
  float theta=0, phi=0, thetai=0, phii=0, theta1=0, phi1=0, theta2=0, phi2=0, beta=0; 
  FILE *fprmap;
  int nmmap=0;
  if (argc >= 5) {

    if ((fprmap = fopen(argv[4], "r")) == NULL) {
        fprintf(stderr, "Error, cannot open input file %s\n", argv[4]);
        return 1;
    } 

    printf("%s loaded!\n", argv[4]);

	while(fgets(line, LINE_LENGTH, fprmap) != NULL){
	    mapid=0; dettype=0; thetai=0; phii=0; theta=0; phi=0; theta1=0; phi1=0; theta2=0; phi2=0; 
        for(i=0; i<LINE_LENGTH; i++){
            if(line[i]=='#'){
                if(PRINT_MAP)printf("%s", line);
                break;
            }
            else if(line[i]>=0){   
	            sscanf(line,"%d\t%c\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", &mapid, &dettype, &detid, &detidi,&beta, &theta, &phi, &thetai, &phii, &theta1, &phi1, &theta2, &phi2);
	            if(PRINT_MAP) printf("%d\t%c\t%d\t%d\t%f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", mapid, dettype, detid, detidi, beta, theta, phi, thetai, phii, theta1, phi1, theta2, phi2);
		        if(mapid >=0 && mapid < MAX_ID) {
                    map2type[mapid] = dettype;
                    map2det[mapid] = detid;
                    map2deti[mapid] = detidi;
		    map2beta[mapid] = beta;
                    mapangles[mapid][0] = theta;
                    mapangles[mapid][1] = phi;
                    mapanglesi[mapid][0]=thetai; 
                    mapanglesi[mapid][1]=phii; 
                    mapangles1[mapid][0]= theta1; 
                    mapangles1[mapid][1]= phi1; 
                    mapangles2[mapid][0]= theta2; 
                    mapangles2[mapid][1]= phi2; 
                    nmmap++;
          		    break;
			    }			    
			    else {
				    printf("Error in reading %s : bad id or format\n", argv[4]);
				    return -1;
			    }
            }
	        else if(line[i]=='\n'){
                if(PRINT_MAP) printf("\n");
                break;
            }
            else {
                continue;
            }

        }
        	memset(line, 0, LINE_LENGTH);
 	  }
 	  
 	fclose(fprmap);
    printf("read %d id maps\n", nmmap);
  }


//GAGG_Calibration_Parameter:


if (argc >=6){
FILE *fgagg;
fgagg=fopen(argv[5],"r");
int idgagg;
int n;
float slopegagg,interceptgagg;

if (fgagg == NULL){
fprintf(stderr,"Error, cannot open input file %s\n",argv[5]);
}
while (fgets(line,LINE_LENGTH,fgagg) != NULL){
sscanf(line,"%d\t%f\t%f\n",&idgagg,&slopegagg,&interceptgagg);
gaid[n]=idgagg;
gaggslope[n]=slopegagg;
gaggintercept[n]=interceptgagg;
n++;
memset(line,0,LINE_LENGTH);
}
fclose(fgagg);
}
  
  
//===========================================================================
//open banana gate file: \\C. W
float banx,bany;
int gaggid,gaggnmpts;
FILE *fprbangate;
char *p;


int lg=0,g1=0,g2=0,g3=0,g4=0,g5=0,g6=0,g7=0,g8=0,g9=0,
    g10=0,g11=0,g12=0,g13=0,g14=0,g15=0,g16=0,g17=0,
    g18=0,g19=0,g20=0,g21=0,
    g22=0,g23=0,g24=0,g25=0,g26=0;

if (argc >= 7){
if ((fprbangate = fopen(argv[6],"r")) == NULL){
        fprintf(stderr, "Error, cannot open input file %s\n", argv[6]);
        return 1;
      }

printf("%s loaded!\n", argv[6]);

while (fgets(line,LINE_LENGTH,fprbangate) !=NULL){
    if ( (p=strchr(line,'#')) != NULL){
      sscanf(line,"%d\t%d",&gaggid,&gaggnmpts);
      GID[lg]=gaggid;
      GNMPTS[lg]=gaggnmpts;
      lg++;
      }
    else{
    sscanf(line,"%d\t%d\t%f\t%f\n",&gaggid,&gaggnmpts,&banx,&bany);
    if (gaggid == GID[0] && gaggnmpts == GNMPTS[0]){
        polyX21[g1]=banx;
        polyY21[g1]=bany;
      printf("%f %f\n",polyX21[g1],polyY21[g1]); 
       g1++;
        }
    if (gaggid == GID[1] && gaggnmpts == GNMPTS[1]){
        polyX22[g2]=banx;
        polyY22[g2]=bany;
      printf("%f %f\n",polyX22[g2],polyY22[g2]); 
       g2++;
         }
    if (gaggid == GID[2] && gaggnmpts == GNMPTS[2]){
        polyX23[g3]=banx;
        polyY23[g3]=bany;
       printf("%f %f\n",polyX23[g3],polyY23[g3]);
        g3++;
          }
    if (gaggid == GID[3] && gaggnmpts == GNMPTS[3]){
        polyX24[g4]=banx;
        polyY24[g4]=bany;
        printf("%f %f\n",polyX24[g4],polyY24[g4]);
        g4++;
          }
    if (gaggid == GID[4] && gaggnmpts == GNMPTS[4]){
        polyX25[g5]=banx;
        polyY25[g5]=bany;
        printf("%f %f\n",polyX25[g5],polyY25[g5]);
        g5++;
          }
    if (gaggid == GID[5] && gaggnmpts == GNMPTS[5]){
        polyX26[g6]=banx;
        polyX26[g6]=bany;
        printf("%f %f\n",polyX26[g6],polyY26[g6]);
        g6++;
          }
    if (gaggid == GID[6] && gaggnmpts == GNMPTS[6]){
        polyX27[g7]=banx;
        polyY27[g7]=bany;
        printf("%f %f\n",polyX27[g7],polyY27[g7]);
        g7++;
          }
    if (gaggid == GID[7] && gaggnmpts == GNMPTS[7]){
        polyX28[g8]=banx;
        polyY28[g8]=bany;
        printf("%f %f\n",polyX28[g8],polyY28[g8]);
        g8++;
          }
    if (gaggid == GID[8] && gaggnmpts == GNMPTS[8]){
        polyX29[g9]=banx;
        polyY29[g9]=bany;
        printf("%f %f\n",polyX29[g9],polyY29[g9]);
        g9++;
          }
    if (gaggid == GID[9] && gaggnmpts == GNMPTS[9]){
        polyX210[g10]=banx;
        polyY210[g10]=bany;
        printf("%f %f\n",polyX210[g10],polyY210[g10]);
        g10++;
          }
    if (gaggid == GID[10] && gaggnmpts == GNMPTS[10]){
        polyX41[g11]=banx;
        polyY42[g11]=bany;
        printf("%f %f\n",polyX41[g11],polyY41[g11]);
        g11++;
          }
    if (gaggid == GID[11] && gaggnmpts == GNMPTS[11]){
        polyX42[g12]=banx;
        polyY42[g12]=bany;
        printf("%f %f\n",polyX42[g12],polyY42[g12]);
        g12++;
          }
    if (gaggid == GID[12] && gaggnmpts == GNMPTS[12]){
        polyX43[g13]=banx;
        polyY43[g13]=bany;
        printf("%f %f\n",polyX43[g13],polyY43[g13]);
        g13++;
          }
    if (gaggid == GID[13] && gaggnmpts == GNMPTS[13]){
        polyX44[g14]=banx;
        polyY44[g14]=bany;
        printf("%f %f\n",polyX44[g14],polyY44[g14]);
        g14++;
          }
    if (gaggid == GID[14] && gaggnmpts == GNMPTS[14]){
        polyX45[g15]=banx;
        polyY45[g15]=bany;
        printf("%f %f\n",polyX45[g15],polyY45[g15]);
        g15++;
          }
    if (gaggid == GID[15] && gaggnmpts == GNMPTS[15]){
        polyX46[g16]=banx;
        polyY46[g16]=bany;
        printf("%f %f\n",polyX46[g16],polyY46[g16]);
        g16++;
          }
    if (gaggid == GID[16] && gaggnmpts == GNMPTS[16]){
        polyX47[g17]=banx;
        polyY47[g17]=bany;
        printf("%f %f\n",polyX47[g17],polyY47[g17]);
        g17++;
          }
    if (gaggid == GID[17] && gaggnmpts == GNMPTS[17]){
        polyX48[g18]=banx;
        polyY48[g18]=bany;
        printf("%f %f\n",polyX48[g18],polyY48[g18]);
        g18++;
          }
    if (gaggid == GID[18] && gaggnmpts == GNMPTS[18]){
        polyX49[g19]=banx;
        polyY49[g19]=bany;
        printf("%f %f\n",polyX49[g19],polyY49[g19]);
        g19++;
          }
    if (gaggid == GID[19] && gaggnmpts == GNMPTS[19]){
        polyX410[g20]=banx;
        polyY410[g20]=bany;
        printf("%f %f\n",polyX410[g20],polyY410[g20]);
        g20++;
          }
    if (gaggid == GID[20] && gaggnmpts == GNMPTS[20]){
        polyX411[g21]=banx;
        polyY411[g21]=bany;
        printf("%f %f\n",polyX411[g21],polyY411[g21]);
        g21++;
          }
    if (gaggid == GID[21] && gaggnmpts == GNMPTS[21]){
        polyX412[g22]=banx;
        polyY412[g22]=bany;
        printf("%f %f\n",polyX412[g22],polyY412[g22]);
        g22++;
          }          
    if (gaggid == GID[22] && gaggnmpts == GNMPTS[22]){
        polyX413[g23]=banx;
        polyY413[g23]=bany;
        printf("%f %f\n",polyX413[g23],polyY413[g23]);
        g23++;
          }
    if (gaggid == GID[23] && gaggnmpts == GNMPTS[23]){
        polyX414[g24]=banx;
        polyY414[g24]=bany;
        printf("%f %f\n",polyX414[g24],polyY414[g24]);
        g24++;
          }
    if (gaggid == GID[24] && gaggnmpts == GNMPTS[24]){
        polyX415[g25]=banx;
        polyY415[g25]=bany;
        printf("%f %f\n",polyX415[g25],polyY415[g25]);
        g25++;
          }
    if (gaggid == GID[25] && gaggnmpts == GNMPTS[25]){
        polyX416[g26]=banx;
        polyY416[g26]=bany;
        printf("%f %f\n",polyX416[g26],polyY416[g26]);
        g26++;
          }
     }
memset(line,0,LINE_LENGTH);
}//end while

fclose(fprbangate);

}//end argc >=7
  
  
  /////////////////////
  // MAIN WHILE LOOP //
  /////////////////////
  while (1) { //main while loop 


      /////////////////////////////////
      // UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////
      etime=-1; tdif=-1; sevtmult=0;  
      //memset(&subevt, 0, sizeof(subevt)); //not needed since everything is redefined (except maybe trace on pileup evts)
      while (1) { //get subevents and event build for one "event" 
        
       // memset(&subevt[sevtmult], 0, sizeof(subevt[sevtmult])); //not needed since everything is redefined (except maybe trace on pileup evts)
        
        //read 4-byte header
        if (fread(sub, sizeof(int)*HEADER_LENGTH, 1, fpr) != 1) break;
        subevt[sevtmult].chn = sub[0] & 0xF;
        subevt[sevtmult].sln = (sub[0] & 0xF0) >> 4;
        subevt[sevtmult].crn = (sub[0] & 0xF00) >> 8;
        subevt[sevtmult].id = subevt[sevtmult].crn*MAX_BOARDS_PER_CRATE*MAX_CHANNELS_PER_BOARD + (subevt[sevtmult].sln-BOARD_START)*MAX_CHANNELS_PER_BOARD + subevt[sevtmult].chn;   
        subevt[sevtmult].hlen = (sub[0] & 0x1F000) >> 12;
        subevt[sevtmult].elen = (sub[0] & 0x7FFE0000) >> 17;
        subevt[sevtmult].fcode = (sub[0] & 0x80000000) >> 31;
        subevt[sevtmult].time = ( (long long int)(sub[2] & 0xFFFF) << 32) + sub[1];
        subevt[sevtmult].ctime = (sub[2] & 0x7FFF0000) >> 16;
        subevt[sevtmult].ctimef = (sub[2] & 0x80000000) >> 31;
        subevt[sevtmult].energy = (sub[3] & 0xFFFF);
        subevt[sevtmult].trlen = (sub[3] & 0x7FFF0000) >> 16;
        subevt[sevtmult].trwlen = subevt[sevtmult].trlen / 2;
        subevt[sevtmult].extra = (sub[3] & 0x80000000) >> 31;       
 
        //rebin raw trap energy from 32k to ....            
        tempf = (float)subevt[sevtmult].energy/RAWE_REBIN_FACTOR;// + RAND;       
        subevt[sevtmult].energy = (int)tempf;                        
 
        //check lengths (sometimes all of the bits for trace length are turned on ...)
       /* if (subevt[sevtmult].elen - subevt[sevtmult].hlen != subevt[sevtmult].trwlen) {
            printf("SEVERE ERROR: event, header, and trace length inconsistencies found\n");
            printf("event length = %d\n", subevt[sevtmult].elen);
            printf("header length = %d\n", subevt[sevtmult].hlen);
            printf("trace length = %d\n", subevt[sevtmult].trwlen);  
            printf("Extra = %d\n", subevt[sevtmult].extra); 
            printf("fcode = %d\n", subevt[sevtmult].fcode);              
            //sleep(1);          
            //return 0;
        } */ 
        
       
        //Set reference time for event building
        if (etime == -1) {
            etime = subevt[sevtmult].time;
            tdif = 0;
        }
        else {
            tdif = subevt[sevtmult].time - etime;
            if (tdif < 0) {
                printf("SEVERE ERROR: tdiff < 0, file must be time sorted\n");
                printf("etime = %lld, time = %lld, and tdif = %lld\n", etime, subevt[sevtmult].time, tdif);                
                return 0;   
            }    
        }    
      
        //Check for end of event, rewind, and break out of while loop
        if (tdif > EVENT_BUILD_TIME) {
            fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR); //fwrite/fread is buffered by system ; storing this in local buffer is no faster!
            break;           
        }    
               
        
        //time between sequential events for a single channel ; useful for determining optimal event building time
        temptime = (subevt[sevtmult].time - idtime[subevt[sevtmult].id])/100; //rebin to 1 micro-second
        if ( temptime >= 0 && temptime < 8192) {
            tdifid[subevt[sevtmult].id][temptime]++;
        }    
        idtime[subevt[sevtmult].id]=subevt[sevtmult].time; //store time for next subevent of channel    
    
        // total pileups
        if (subevt[sevtmult].fcode==1) {
            pileupcount++;
        }
        
        //Histogram raw spectra
        hit[0][subevt[sevtmult].id]++;
        if (subevt[sevtmult].fcode==1) 
        hit[1][subevt[sevtmult].id]++;     
  
        if (subevt[sevtmult].energy >= 0 && subevt[sevtmult].energy < 8192)
        e_raw[subevt[sevtmult].id][subevt[sevtmult].energy]++;
 
        if (subevt[sevtmult].time/1000000000 >= 0 && subevt[sevtmult].time/1000000000 < 8192) // rebin to 10 seconds
        tevt_raw[subevt[sevtmult].id][subevt[sevtmult].time/1000000000]++; // rebin to 10 seconds

        if (subevt[sevtmult].ctime >= 0 && subevt[sevtmult].ctime < 8192)
        tcfd_raw[subevt[sevtmult].id][subevt[sevtmult].ctime]++;

        if (tdif >= 0 && tdif < 4096 && sevtmult!=0)
        tdif_raw[subevt[sevtmult].id][tdif]++;


        //if CFD is enabled, ctime will be non-zero
        //tempf = (float)subevt[sevtmult].ctime*10.0/32768.0;
        //subevt[sevtmult].time = subevt[sevtmult].time + (long long int)tempf;    

        //Calibrate energy and time
        tempf = ((float)subevt[sevtmult].energy*ecal[subevt[sevtmult].id][1] + ecal[subevt[sevtmult].id][0])/new_gain;// + RAND;       
        subevt[sevtmult].energy = (int)tempf; 
        //subevt[sevtmult].time += (long long int)tcal[subevt[sevtmult].id][0];       
	
        //Histogram calibrated spectra
        if (subevt[sevtmult].energy >= 0 && subevt[sevtmult].energy < 8192)
        e_cal[subevt[sevtmult].id][subevt[sevtmult].energy]++;        
        
        if (subevt[sevtmult].time/1000000000 >= 0 && subevt[sevtmult].time/1000000000 < 8192)        
        tevt_cal[subevt[sevtmult].id][subevt[sevtmult].time/1000000000]++;
         
        //continue on if no trace, esum, or qsum
        if (subevt[sevtmult].hlen==HEADER_LENGTH && subevt[sevtmult].trwlen==0 ) {
            sevtmult++;
            continue;
        }
        
        //more data than just the header; read entire sub event
        fseek(fpr, -sizeof(int)*HEADER_LENGTH, SEEK_CUR);
        if (fread(sub, sizeof(int)*subevt[sevtmult].elen, 1, fpr) != 1) break;
                              
        //trace
        k=0;
        for (i = subevt[sevtmult].hlen; i < subevt[sevtmult].elen; i++) {      
            subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0x3FFF; // the upper 2 bits/16 bits are filled with 0s
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k] = sub[i] & 0xFFFF; //C. W  
           subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0x3FFF;
           // subevt[sevtmult].tr[i - subevt[sevtmult].hlen + k + 1] = (sub[i]>>16) & 0xFFFF; //C. W
            k=k+1;
        } 
        
     // if (subevt[sevtmult].id == 4 && subevt[sevtmult].fcode == 1) DB(subevt[sevtmult].tr);
            
        //continue if no esum or qsum   
        if (subevt[sevtmult].hlen==HEADER_LENGTH) {
            sevtmult++;        
            continue;
        }
        
        //esum
        if (subevt[sevtmult].hlen==8 || subevt[sevtmult].hlen==16) { 
            for (i=4; i < 8; i++) {
                subevt[sevtmult].esum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[sevtmult].hlen==12) { 
            for (i=4; i < 12; i++) {
                subevt[sevtmult].qsum[i-4] = sub[i];
            }
        }
    
        //qsum
        if (subevt[sevtmult].hlen==16) { 
            for (i=8; i < 16; i++) {
                subevt[sevtmult].qsum[i-8] = sub[i];
            }
        }    
        
        sevtmult++;
     
      } //end while loop for unpacking sub events and event building for one "event"
      if (sevtmult==0) break; //end main WHILE LOOP when out of events 
      mult[0][sevtmult]++; //Histogram raw sub event multiplicity 
      sevtcount += sevtmult;
      evtcount++; //event-built number
      /////////////////////////////////////
      // END UNPACK DATA AND EVENT BUILD //
      /////////////////////////////////////


      //skip detector building below if no map file
      if (argc >= 5) {        


      ////////////////////////////////////////
      // MAP SUB EVENTS INTO DETECTOR TYPES //
      ////////////////////////////////////////
      memset(&ge, 0, sizeof(ge)); //This is needed but could be replaced by setting suppress, pileup, nonprompt, clean, x/s/bmult to zero at start of loop!
      memset(&si,0,sizeof(si)); //C. W
    //Start Mapping Raw Subevents into Detector Types:
    detmaps(sevtmult, subevt, ge, si);
      ////////////////////////////////////////////
      // END MAP SUB EVENTS INTO DETECTOR TYPES //
      ////////////////////////////////////////////

      ////////////////////////////
      // PROCESS DETECTOR TYPES //
      ////////////////////////////
//C.W
//printf("hello\n");

//Only Generate PID when Germanium Energy conditions meet:
/*
int km;
int gevalid=0;
//For each subevent search for the gamma ray energies that meet the gating condition to evaluate whether there is a corresponding particle and gamma.
for (km=0;km<gmult;km++){
  if ((ge[km].energy > Emin && ge[km].energy < Emax)) {
  gevalid++;
  }
}
*/

/////////////////////////
// Si Detector Type /////
/////////////////////////    //C. W

//Checking Gagg and Germanium Multiplicities:
int gaggvalid=gaggproc(si);
gmult=gamproc(ge, ge_spe_xtl, ge_spe);

//Checking the energy of GAGG thas fire more than twices for a given events:

int p;
if (gaggvalid > 0  && gmult>0){
  for (p=0;p<gaggvalid;p++){
    if ((si[p].energy < 200)) 
    printf("gam: %d,gem: %d,gaggenergy: %d, tpr: %d, id: %d\n",gaggvalid,gmult,(int)(si[p].energy/100.),(int)(si[p].tpratio),si[p].id);
                  }
}
//=====================

//Perform Kinematic Corrections:
doppcorrect(gaggvalid, gmult, si, ge, option, ml, mbeam, mh, vbz, betascale);


//Record Particle Hits:
parthit(gaggvalid, gmult, ge, si, e_lit2, e_lit4, theta_res2, theta_res4, gagghist);


//Perform 1D gagg-Gamma Histogram:
gammagagg1Dspectra(gaggvalid, gmult, si, ge, gagg_gamma_tdiff_ring2, gagg_gamma_tdiff_ring4, gagg_gamma_tdiff2,gated_gamma, gated_gammatrue,  gated_gamma_np, gated_gamma_nptrue);



/*

//C.W
//if (bantest(si[1].tail,si[1].peak,polyX,polyY,24) == 1 || bantest(si[2].tail,si[2].peak,polyX2,polyY2,28) == 1
//	|| bantest(si[3].tail,si[3].peak,polyX3,polyY3,26) == 1 || bantest(si[4].tail,si[4].peak,polyX4,polyY4,24) ==1){
 

*/

//} //end if condition for GAGG+Gamma Gamma

///////////////////////////////////////////////////////////C.W
//1D Histogram Rough Gate:
/*
int gemdet;
if (gaggvalid > 0 && gmult > 0){
    for (gemdet=0;gemdet<gmult;gemdet++){
        if (ge[gemdet].energy > 0 && ge[gemdet].energy < 4096)
            ge_gagg_spe_r[ge[gemdet].id][ge[gemdet].energy]++;
                                        }
                               }
*/
////////////////////////////////////////////////////////////////

      ///////////////////////////
      // END PROCESS DETECTORS //
      ///////////////////////////



      ////////////////////////
      // FINAL USER SPECTRA //
      ////////////////////////

      
//if (bantest(si[1].tail,si[1].peak,polyX,polyY,24) == 1 || bantest(si[2].tail,si[2].peak,polyX2,polyY2,28) == 1
 // || bantest(si[3].tail,si[3].peak,polyX3,polyY3,26) == 1 || bantest(si[4].tail,si[4].peak,polyX4,polyY4,24) ==1)


//Only gg that come coincidences with proton allowed
/*
if(gaggvalid == 1 && gmult > 0){
long long int tdif1,tdif2;
      //gamma-gamma time and energy 
      for (i=0; i<gmult; i++) {
        for (j=0; j<gmult; j++) {
            if (i!=j) {
                if (ge[i].validp == 1 && ge[j].validp == 1)
                    tdif1 = ge[i].time - ge[j].time + 2000;
                if (ge[i].validnp == 1 && ge[j].validnp == 1)
                    tdif2 = ge[i].time - ge[j].time + 2000; 
                //time difference matrix between ge-ge within prompt region:    
                if (tdif1 >= 0 && tdif1 < 4096 && ge[i].energy >= 0 && ge[i].energy < 4096 && ge[j].energy >= 0 && ge[j].energy < 4096 && (ge[i].validp == 1) && ge[j].validp == 1) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif1[ge[i].energy][tdif1]++;
                    else 
                    gg_tdif1[ge[j].energy][tdif1]++;
                }   
               //time difference matrix between ge-ge within non-prompt region:
                if (tdif2 >=0 && tdif2 < 4096 && ge[i].energy >=0 && ge[j].energy < 4096 && ge[j].energy >=0 && ge[j].energy < 4096 && ge[i].validnp == 1 && ge[j].validnp == 1) {
                    if (ge[i].energy < ge[j].energy)
                    gg_tdif2[ge[i].energy][tdif2]++;
                    else
                    gg_tdif2[ge[j].energy][tdif2]++;
                }


                //prompt gamma-gamma
               // if (tdif1 >= 1967 && tdif1 <= 2032 && ge[i].validp == 1 && ge[j].validp == 1) {
                if (ge[i].validp == 1 && ge[j].validp == 1) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_prompt[ge[i].energy][ge[j].energy]++;
		          if(ge[i].etrue >=0 && ge[i].etrue < 5000 && ge[j].energy >=0 && ge[j].energy < 5000){
               gg_prompt_nddopp[(int)ge[i].etrue][ge[j].energy]++;} //C.W
              if (ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){
               gg_prompt_ddopp[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W
                }        
                
                //prompt-non prompt:
                
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validp == 1 && ge[j].validnp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt1[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp1[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
              if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp1[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validnp == 1 && ge[j].validp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt2[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp2[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp2[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                
                
                
                //non-prompt and prompt:
                 //non-prompt gamma-gamma (mult by 0.2555)
               // if ( (tdif2 >= 1963 && tdif2 <= 2040) && ge[i].validnp == 1 && ge[j].validnp == 1 ) {
                if ( ge[i].validnp == 1 && ge[j].validnp == 1 ) {

                    if (ge[i].energy >= 0 && ge[i].energy < 5000 && ge[j].energy >= 0 && ge[j].energy < 5000) 
                    gg_nprompt[ge[i].energy][ge[j].energy]++;
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].energy >=0 && ge[j].energy <5000){
            		    gg_nprompt_nddopp[(int)ge[i].etrue][ge[j].energy]++;} //C.W 
                if(ge[i].etrue >=0 && ge[i].etrue <5000 && ge[j].etrue >=0 && ge[j].etrue <5000){ 
                   gg_nprompt_ddopp[(int)ge[i].etrue][(int)ge[j].etrue]++;} //C.W              
                }    
                           
            }
        }
      }  
}//end gagg_valid

*/



pol(gaggvalid, gmult,ge,pol_para, pol_perp);



      ////////////////////////////
      // END FINAL USER SPECTRA //
      ////////////////////////////

      } //end argc >= 5 condition 

      //event stats, print status every 10000 events
      lle_div=lldiv(evtcount,10000);
      if ( lle_div.rem == 0 && strstr(argv[1], "p") != NULL) {
        fprpos = ftell(fpr);
        tempf = (float)fprsize/(1024.*1024.*1024.);
        printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\r", sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
      }      

  
//printf("Ep<200: %d,Ep>200: %d\n",plt200,pgt200);
        
  } // end main while loop 
  /////////////////////////
  // END MAIN WHILE LOOP //
  /////////////////////////
  fprpos = ftell(fpr);
  tempf = (float)fprsize/(1024.*1024.*1024.);
  printf("Total SubEvents: \x1B[32m%llu \x1B[31m(%d%% pileup)\x1B[0m\nTotal Events: \x1B[32m%llu (%.1f <mult>)\x1B[0m\nPercent Complete: \x1B[32m%ld%% of %.3f GB\x1B[0m\n\033[3A\r", sevtcount, (int)((100*pileupcount)/sevtcount), evtcount, (float)sevtcount/(float)evtcount, (100*fprpos/fprsize), tempf);
           

  
  
  
  
  
  ////////////////////
  // WRITE SPECTRA //
  ///////////////////
  printf("\n\n\n\nWriting Spectra to Disk ...\n");
     
  //Event Spectra   
  /*
  write_data4(HIT, *hit, 4096, 2, overwrite);
  write_data4(MULT, *mult, 4096, 1, overwrite);
  write_data4(TDIFID, *tdifid, MAX_ID, 8192, overwrite);
  write_data4(E_RAW, *e_raw, MAX_ID, 8192, overwrite);
  write_data4(E_CAL, *e_cal, MAX_ID, 8192, overwrite);
  write_data4(TEVT_RAW, *tevt_raw, MAX_ID, 8192, overwrite);
  write_data4(TEVT_CAL, *tevt_cal, MAX_ID, 8192, overwrite);
  write_data4(TCFD_RAW, *tcfd_raw, MAX_ID, 8192, overwrite);
  write_data4(TDIF_RAW, *tdif_raw, MAX_ID, 4096, overwrite);
  write_data4(TDIF_CAL, *tdif_cal, MAX_ID, 4096, overwrite);
  write_data4(TDIF_CAL0_ETHRESH, *tdif_cal0_ethresh, MAX_ID, 4096, overwrite);
  write_data4(IDID, *idid, MAX_ID, MAX_ID, overwrite);

*/
  //Detector Processed Spectra
  //Ge
 /*
  write_data4(GE_BGO_TDIF, *ge_bgo_tdif, MAX_GE, 4096, overwrite);
  write_data4(GE_XTL_TDIF, *ge_xtl_tdif, MAX_GE, 4096, overwrite);
  write_data4(GE_XTL_TDIF_ETHRESH, *ge_xtl_tdif_ethresh, MAX_GE, 4096, overwrite);
  write_data4(GE_SPE_XTL, *ge_spe_xtl, MAX_GE*MAX_GE_XTL, 8192, overwrite);
  write_data4(GE_SPE, *ge_spe, MAX_GE, 8192, overwrite);
  write_data4(GE_SPE_DOPP, *ge_spe_dopp, MAX_GE, 8192, overwrite); //C.W
  write_data4(GE_SPE_CLEAN, *ge_spe_clean, MAX_GE, 8192, overwrite);
  //trinity
  write_data4(PID, *pid, 4096, 4096, overwrite);
  write_data4(PID_EVSP, *pid_evsp, 4096, 4096, overwrite);
  write_data4(PID_EVST, *pid_evst, 4096, 4096, overwrite);
  write_data4(PID_EVSTAU, *pid_evstau, 4096, 4096, overwrite);
  write_data4(PID_EVSR, *pid_evsr, 4096, 4096, overwrite);
  write_data4(PID_TAUVSR, *pid_tauvsr, 4096, 4096, overwrite);
  write_data4(GAGGDT, *gaggdt, 4096, 4096, overwrite);
*/

     write_data4(GE_SPE_CLEAN, *ge_spe_clean, MAX_GE, 4096, overwrite);  
    write_data4(PID_Gamma_R, *ge_gagg_spe_r, MAX_GE, 4096, overwrite);
    write_data4(PID_Gammaa, *gagg_gamma_tdiff_ring2, 4096, 4096, overwrite);
    write_data4(PID_Gammab, *gagg_gamma_tdiff_ring4, 4096, 4096, overwrite);
    write_data4(PID_Gamma2, *gagg_gamma_tdiff2, 5000, 5000, overwrite);
    write_data4(PID_Gamma_1Dp, *gated_gamma, 1, 5000, overwrite);
    write_data4(PID_Gamma_1Dptrue, *gated_gammatrue,1,5000,overwrite);
    write_data4(PID_Gamma_1Dnp, *gated_gamma_np, 1, 5000, overwrite);
    write_data4(PID_Gamma_1Dnptrue, *gated_gamma_nptrue, 1, 5000, overwrite);

/*
  for (int i=0; i<26; ++i) {
    char name[128];
    sprintf(name, "gagg_pid_%i.spn", i+1);
    write_data4(name, &pid_all[i][0][0], 4096, 4096, overwrite);
  }

*/
 
 /*
 write_data4(PID_QDC21, *pid_qdc21, 4096,4096, overwrite);
 write_data4(PID_QDC22, *pid_qdc22, 4096, 4096, overwrite);
 write_data4(PID_QDC23, *pid_qdc23, 4096, 4096, overwrite);
 write_data4(PID_QDC24,*pid_qdc24, 4096, 4096, overwrite);
 write_data4(PID_QDC25, *pid_qdc25, 4096,4096, overwrite); 
 write_data4(PID_QDC26, *pid_qdc26, 4096, 4096, overwrite);
 write_data4(PID_QDC27, *pid_qdc27, 4096, 4096, overwrite);
 write_data4(PID_QDC28,*pid_qdc28, 4096, 4096, overwrite);
 write_data4(PID_QDC29, *pid_qdc29, 4096,4096, overwrite);
 write_data4(PID_QDC210, *pid_qdc210, 4096, 4096, overwrite);
 write_data4(PID_QDC41, *pid_qdc41, 4096, 4096, overwrite);
 write_data4(PID_QDC42,*pid_qdc42, 4096, 4096, overwrite);
 write_data4(PID_QDC43, *pid_qdc43, 4096,4096, overwrite);
 write_data4(PID_QDC44, *pid_qdc44, 4096, 4096, overwrite);
 write_data4(PID_QDC45, *pid_qdc45, 4096, 4096, overwrite);
 write_data4(PID_QDC46,*pid_qdc46, 4096, 4096, overwrite);
 write_data4(PID_QDC47, *pid_qdc47, 4096, 4096, overwrite);
 write_data4(PID_QDC48,*pid_qdc48, 4096, 4096, overwrite);
 write_data4(PID_QDC49, *pid_qdc49, 4096,4096, overwrite);
 write_data4(PID_QDC410, *pid_qdc410, 4096, 4096, overwrite);
 write_data4(PID_QDC411, *pid_qdc411, 4096, 4096, overwrite);
 write_data4(PID_QDC412,*pid_qdc412, 4096, 4096, overwrite);
 write_data4(PID_QDC413, *pid_qdc413, 4096, 4096, overwrite);
 write_data4(PID_QDC414, *pid_qdc414, 4096, 4096, overwrite);
 write_data4(PID_QDC415, *pid_qdc415, 4096, 4096, overwrite);
 write_data4(PID_QDC416, *pid_qdc416, 4096, 4096, overwrite);
 write_data4(GAGG_GE_TDIFF, *gagg_ge_tdiff, 4096, 4096, overwrite);
*/
//write_data4(PID_TRACEINT, *pid_traceint, 4096, 4096, overwrite);
//write_data4(PID_E1, *pid_e1,4096,4096, overwrite);
//write_data4(PID_E2, *pid_e2,4096,4096, overwrite);
/*
  for (i=0; i<500; i++) {
  	if (strace[i]>0) sptrace[0][i] = strace[i]/dbcount;
  }
  write_data4(SPTRACE, *sptrace, 1, 4096, overwrite);

*/
  //Final User Spectra
  //gamma-gamma
  write_data4(GAGG_HIST, *gagghist, 4096, 1, overwrite);
  write_data4(THETA_RES2, *theta_res2, 4096, 28, overwrite);
  write_data4(THETA_RES4, *theta_res4, 4096, 28, overwrite);
  write_data4(E_LIT2, *e_lit2, 4096, 28, overwrite);
  write_data4(E_LIT4, *e_lit4, 4096, 28, overwrite);
  write_data4(POL_PARA,*pol_para, 5000, 5000, overwrite);
  write_data4(POL_PERP, *pol_perp, 5000, 5000, overwrite);
/*
  write_data4(GG_TDIF, *gg_tdif, 4096, 4096, overwrite);
  write_data4(GG_TDIF1, *gg_tdif1, 4096, 4096, overwrite);
  write_data4(GG_TDIF2, *gg_tdif2, 4096, 4096, overwrite);
  write_data4(GG_PROMPT, *gg_prompt, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PURE, *gg_prompt_pure, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PURE, *gg_nprompt_pure, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT, *gg_nprompt, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT1, *gg_nprompt1, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT2, *gg_nprompt2, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_nddopp, *gg_prompt_nddopp, 5000, 5000, overwrite); //C.W
  write_data4(GG_PROMPT_ddopp, *gg_prompt_ddopp, 5000, 5000, overwrite); //C. W
  write_data4(GG_NPROMPT_nddopp, *gg_nprompt_nddopp, 5000,5000, overwrite); //C.W
  write_data4(GG_NPROMPT_nddopp1, *gg_nprompt_nddopp1, 5000,5000, overwrite); //C.W
  write_data4(GG_NPROMPT_ddopp1, *gg_nprompt_ddopp1, 5000, 5000, overwrite); //C.W
  write_data4(GG_NPROMPT_nddopp2, *gg_nprompt_nddopp2, 5000,5000, overwrite); //C.W
  write_data4(GG_NPROMPT_ddopp2, *gg_nprompt_ddopp2, 5000, 5000, overwrite); //C.W
  write_data4(GG_NPROMPT_ddopp, *gg_nprompt_ddopp, 5000, 5000, overwrite);//C.W
 */
  fclose(fpr);
  fclose(debugfile);
  
  return 0;
}


// C function for ban test
//data x   data y   ban x array   ban y array   num ban points

bool bantest(float x, float y, float *polyX, float *polyY, int polyCorners) {

   int   i, j=polyCorners-1 ;
   bool  oddNodes=FALSE      ;

   for (i=0; i<polyCorners; i++) {
     if (( (polyY[i]< y && polyY[j]>=y) || (polyY[j]< y && polyY[i]>=y) ) && (polyX[i]<=x || polyX[j]<=x)) {
         oddNodes^=(polyX[i]+(y-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<x); 
     }    
     j=i; 
   }

  return oddNodes; 
 }
 
 
 
///////////////////////
// Write 2-byte data //
///////////////////////
void write_data2(char *filename, short *data, int xdim, int ydim, int overwrite) { //2byte per channel Write / Add to previous
    
    FILE *FP;
    int i;
    short *previous;
    
    if(!overwrite) {
        //allocate memory for 1d-array for reading in rows of 2d Radware matrix
        if ( ( previous = (short *)malloc(xdim * ydim * sizeof(short)) ) == NULL ) {
            printf("\nError, memory not allocated.\n");
            exit(1);
        }
      
        //open previous spectra file  
        if( (FP=fopen(filename, "r")) != NULL ){
            fread(previous, sizeof(short)*xdim*ydim, 1, FP);        
            fclose(FP);
            //update spectra
            for (i=0; i<xdim*ydim; i++) {
                if(previous[i] < (powf(2,sizeof(short)*8.0)-2)) 
                    data[i] = data[i] + previous[i];
                }   
            }
        else{
            printf("%s did not previously exist, creating ...\n", filename);       
        }     
   
        //Deallocate previous data
        free(previous);
    }
  
    FP=fopen(filename, "w");
    fwrite(data, sizeof(short)*xdim, ydim, FP);
    fclose(FP);            
}    


///////////////////////
// Write 4-byte data //
///////////////////////
void write_data4(char *filename, int *data, int xdim, int ydim, int overwrite) { //4byte per channel Write / Add to previous
    
    FILE *FP;
    int i;
    int *previous;
   
    if(!overwrite) {
        //allocate memory for 1d-array for reading in rows of 2d Radware matrix
        if ( ( previous = (int *)malloc(xdim * ydim * sizeof(int)) ) == NULL ) {
            printf("\nError, memory not allocated.\n");
            exit(1);
        }
      
        //open previous spectra file  
        if( (FP=fopen(filename, "r")) != NULL ){
            fread(previous, sizeof(int)*xdim*ydim, 1, FP);        
            fclose(FP);
            //update spectra
            for (i=0; i<xdim*ydim; i++) {
                if(previous[i] < (powf(2,sizeof(int)*8.0)-2)) 
                    data[i] = data[i] + previous[i];
            }   
        } 
        else{
            printf("%s did not previously exist, creating ...\n", filename);       
        }   
       
        //Deallocate previous data
        free(previous);
    }
  
    FP=fopen(filename, "w");
    fwrite(data, sizeof(int)*xdim, ydim, FP);
    fclose(FP);            
}

void detmaps(int sevtmult, struct subevent *subevt, struct gdetector *ge, struct sidetector *si){
int i,j,l;     
 for (i=0; i<sevtmult; i++) { 
        
/* //07/31/'22 C.W

        for (j=0; j<sevtmult; j++) {
            if (i!=j)
            idid[subevt[i].id][subevt[j].id]++;
        }

*/
        
        //printf("i=%d, sevtmult=%d, subevt[i].id=%d, map2type[subevt[i].id]=%c, subevt[i].energy=%d, subevt[i].time=%lld\n", i, sevtmult, subevt[i].id, map2type[subevt[i].id], subevt[i].energy, subevt[i].time);
        //fflush(stdout);
          
        //if CFD is enabled, ctime will be non-zero
        //tempf = (float)subevt[i].ctime*10.0/32768.0;
        //tempf = (float)subevt[i].ctime*1.0/32768.0; //should make no difference ; need decimals or convert time to 1 ns bins
        //printf("%lld e=%d ", subevt[i].time, subevt[i].energy);
        //subevt[i].time = subevt[i].time + (long long int)tempf;  
        //printf("%lld and %f\n", subevt[i].time, tempf); 
  
               
        //Histogram calibrated tdif spectra (Do here to keep subevents within an event time ordered during event build above)
       /* 07/31/'22 //C.W
	 if (i==0) etime = subevt[i].time;
        tdif = abs(subevt[i].time - etime); 
        if (tdif >= 0 && tdif < 4096 && i!=0) tdif_cal[subevt[i].id][tdif]++;
      */
	  //if (tdif > 20 && tdif < 50) DB(subevt[i].tr);
        
        //tdif with respect to channel id 0
     /*  07/31/'22 C.2
	   if (i!=0 && subevt[0].id==0 && subevt[0].energy >= 10 && subevt[0].energy <= 8000 && subevt[i].energy >= 10  && subevt[i].energy <= 8000) {
            if (tdif >= 0 && tdif < 4096 ) tdif_cal0_ethresh[subevt[i].id][tdif]++;   
        }            
      */  

        /////////////////////
        // G = Ge Detector //
        /////////////////////              
        if ( map2type[subevt[i].id] == 'G' ) { //Keep G and ge or switch to C and clover/cl?
        
            if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_GE_XTL) { //Ge crystal 
               if ( ge[map2det[subevt[i].id]].xmult >= MAX_GE_XTL) {               
                    printf("SEVERE ERROR: Same Ge(xtl) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].xid[ge[map2det[subevt[i].id]].xmult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].xe[ge[map2det[subevt[i].id]].xmult] = subevt[i].energy; 
		ge[map2det[subevt[i].id]].xedopp[ge[map2det[subevt[i].id]].xmult] = subevt[i].energy/(1+map2beta[subevt[i].id]*cos(PIR*mapanglesi[subevt[i].id][0])); //C .W              
                ge[map2det[subevt[i].id]].xt[ge[map2det[subevt[i].id]].xmult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].xct[ge[map2det[subevt[i].id]].xmult] = subevt[i].ctime;  
                ge[map2det[subevt[i].id]].xpileup[ge[map2det[subevt[i].id]].xmult] = subevt[i].fcode;  
                ge[map2det[subevt[i].id]].xsubevtid[ge[map2det[subevt[i].id]].xmult] = i;  
                        
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][0] = mapangles[subevt[i].id][0];   
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][1] = mapanglesi[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][2] = mapangles1[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xtheta[ge[map2det[subevt[i].id]].xmult][3] = mapangles2[subevt[i].id][0];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][0] = mapangles[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][1] = mapanglesi[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][2] = mapangles1[subevt[i].id][1];               
                ge[map2det[subevt[i].id]].xphi[ge[map2det[subevt[i].id]].xmult][3] = mapangles2[subevt[i].id][1]; 

                ge[map2det[subevt[i].id]].xmult++;   
            }                                                                
            if (map2deti[subevt[i].id] > MAX_GE_XTL && map2deti[subevt[i].id] <= MAX_GE_XTL + MAX_GE_SEG) { //Ge segment
               if ( ge[map2det[subevt[i].id]].smult >= MAX_GE_SEG ) {               
                    printf("SEVERE ERROR: Same Ge(seg) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].sid[ge[map2det[subevt[i].id]].smult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].se[ge[map2det[subevt[i].id]].smult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].st[ge[map2det[subevt[i].id]].smult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].sct[ge[map2det[subevt[i].id]].smult] = subevt[i].ctime;  
                ge[map2det[subevt[i].id]].spileup[ge[map2det[subevt[i].id]].smult] = subevt[i].fcode;  
                ge[map2det[subevt[i].id]].ssubevtid[ge[map2det[subevt[i].id]].smult] = i;  
                            
                ge[map2det[subevt[i].id]].smult++;               
            }            
            if (map2deti[subevt[i].id] > MAX_GE_XTL + MAX_GE_SEG) { //BGO 
               if ( ge[map2det[subevt[i].id]].bgomult >= MAX_GE_BGO ) {               
                    printf("SEVERE ERROR: Same Ge(bgo) twice within event build; Make event build time smaller!\n");    
                    printf("ge[map2det[subevt[i].id]].xmult=%d, ge[map2det[subevt[i].id]].smult=%d, ge[map2det[subevt[i].id]].bgomult=%d\n", ge[map2det[subevt[i].id]].xmult, ge[map2det[subevt[i].id]].smult,ge[map2det[subevt[i].id]].bgomult);
                    continue;
                    //return -1;
                }   
                ge[map2det[subevt[i].id]].bgoid[ge[map2det[subevt[i].id]].bgomult] = map2deti[subevt[i].id];               
                ge[map2det[subevt[i].id]].bgoe[ge[map2det[subevt[i].id]].bgomult] = subevt[i].energy;               
                ge[map2det[subevt[i].id]].bgot[ge[map2det[subevt[i].id]].bgomult] = subevt[i].time;                           
                ge[map2det[subevt[i].id]].bgoct[ge[map2det[subevt[i].id]].bgomult] = subevt[i].ctime;              
                ge[map2det[subevt[i].id]].bgopileup[ge[map2det[subevt[i].id]].bgomult] = subevt[i].fcode;
                ge[map2det[subevt[i].id]].bgosubevtid[ge[map2det[subevt[i].id]].bgomult] = i;  

                ge[map2det[subevt[i].id]].bgomult++;               
            } 
                 
                
        } //end G

        ///////////////////           //C. W
        // S=Si Detector //
	if ( map2type[subevt[i].id] == 'S' && subevt[i].trlen == 0 ){
	if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_SI_PAIR){
		            si[map2det[subevt[i].id]].siid[si[map2det[subevt[i].id]].simult] = map2deti[subevt[i].id];               
                si[map2det[subevt[i].id]].sie[si[map2det[subevt[i].id]].simult] = subevt[i].energy;             
                si[map2det[subevt[i].id]].sit[si[map2det[subevt[i].id]].simult] = subevt[i].time;                           
                si[map2det[subevt[i].id]].sict[si[map2det[subevt[i].id]].simult] = subevt[i].ctime;  
                si[map2det[subevt[i].id]].sipileup[si[map2det[subevt[i].id]].simult] = subevt[i].fcode;  
                si[map2det[subevt[i].id]].sisubevtid[si[map2det[subevt[i].id]].simult] = i;        
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][0];   
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][0];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][1]; 
		            si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] = subevt[i].qsum[0];   
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] = subevt[i].qsum[1];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] = subevt[i].qsum[2];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] = subevt[i].qsum[3];  
		            si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] = subevt[i].qsum[4];   
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] = subevt[i].qsum[5];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] = subevt[i].qsum[6];               
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] = subevt[i].qsum[7];              
                si[map2det[subevt[i].id]].simult++;   
                
		}

} //end Si        //C. W

// Another types of GAGG-mappings:
if ( map2type[subevt[i].id] == 'S' && subevt[i].trlen != 0 ){
	if (map2deti[subevt[i].id] > 0 && map2deti[subevt[i].id] <= MAX_SI_PAIR){
		            si[map2det[subevt[i].id]].siid[si[map2det[subevt[i].id]].simult] = map2deti[subevt[i].id];               
                si[map2det[subevt[i].id]].sie[si[map2det[subevt[i].id]].simult] = subevt[i].energy;             
                si[map2det[subevt[i].id]].sit[si[map2det[subevt[i].id]].simult] = subevt[i].time;                           
                si[map2det[subevt[i].id]].sict[si[map2det[subevt[i].id]].simult] = subevt[i].ctime;  
                si[map2det[subevt[i].id]].sipileup[si[map2det[subevt[i].id]].simult] = subevt[i].fcode;  
                si[map2det[subevt[i].id]].sisubevtid[si[map2det[subevt[i].id]].simult] = i;        
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][0];   
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][0];               
                si[map2det[subevt[i].id]].sitheta[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][0];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][0] = mapangles[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][1] = mapanglesi[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][2] = mapangles1[subevt[i].id][1];               
                si[map2det[subevt[i].id]].siphi[si[map2det[subevt[i].id]].simult][3] = mapangles2[subevt[i].id][1]; 
//Reconstructing Trace to Qsum: //C.W
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] = 0;
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] = 0;

for (l=0; l<31; l++){
                si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][0] + subevt[i].tr[l];
}   

for (l=31;l<60;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][1] + subevt[i].tr[l];                    
}

for (l=60;l<75;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][2] + subevt[i].tr[l];               
}

for (l=75;l<95;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][3] + subevt[i].tr[l];
}  

for (l=95;l<105;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][4] + subevt[i].tr[l];   
}

for (l=105;l<160;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][5] + subevt[i].tr[l];
}

for (l=160;l<175;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][6] + subevt[i].tr[l];               
}

for (l=175;l<200;l++){
si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] = si[map2det[subevt[i].id]].siqdc[si[map2det[subevt[i].id]].simult][7] + subevt[i].tr[l];              
}

si[map2det[subevt[i].id]].simult++;   
                
		}

} //end Si        //C. W


        ///////////////////  
        
      } // end i loop over sevtmult

}     

int gaggproc(struct sidetector *si){

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


si[i].traceint=(((si[i].siqdc[0][2]+si[i].siqdc[0][3]+si[i].siqdc[0][4]+si[i].siqdc[0][5]+si[i].siqdc[0][6]-(115./(31.+29.))*(si[i].siqdc[0][0]+si[i].siqdc[0][1]))+(si[i].siqdc[1][2]+si[i].siqdc[1][3]+si[i].siqdc[1][4]+si[i].siqdc[1][5]+si[i].siqdc[1][6]-(115./(31.+29.))*(si[i].siqdc[1][0]+si[i].siqdc[1][1])))/30.);

si[i].tpratio=4000.*si[i].tail/si[i].peak;

si[i].id = i;
si[i].valid = 1;

if (si[i].id == gaid[i-1]){
si[i].energy=100.*(gaggslope[i-1]*si[i].traceint+gaggintercept[i-1]);
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
    if (i == 1) { ingate = bantest(si[1].tpratio,si[1].energy,polyX21,polyY21,GNMPTS[0]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 2) { ingate = bantest(si[2].tpratio,si[2].energy,polyX22,polyY22,GNMPTS[1]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 3) { ingate = bantest(si[3].tpratio,si[3].energy,polyX23,polyY23,GNMPTS[2]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 4) { ingate = bantest(si[4].tpratio,si[4].energy,polyX24,polyY24,GNMPTS[3]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 5) { ingate = bantest(si[5].tpratio,si[5].energy,polyX25,polyY25,GNMPTS[4]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 6) { ingate = bantest(si[6].tpratio,si[6].energy,polyX26,polyY26,GNMPTS[5]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 7) { ingate = bantest(si[7].tpratio,si[7].energy,polyX27,polyY27,GNMPTS[6]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 8) { ingate = bantest(si[8].tpratio,si[8].energy,polyX28,polyY28,GNMPTS[7]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 9) { ingate = bantest(si[9].tpratio,si[9].energy,polyX29,polyY29,GNMPTS[8]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 10) { ingate = bantest(si[10].tpratio,si[10].energy,polyX210,polyY210,GNMPTS[9]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 11) { ingate = bantest(si[11].tpratio,si[11].energy,polyX41,polyY41,GNMPTS[10]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 12) { ingate = bantest(si[12].tpratio,si[12].energy,polyX42,polyY42,GNMPTS[11]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 13) { ingate = bantest(si[13].tpratio,si[13].energy,polyX43,polyY43,GNMPTS[12]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 14) { ingate = bantest(si[14].tpratio,si[14].energy,polyX44,polyY44,GNMPTS[13]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 15) { ingate = bantest(si[15].tpratio,si[15].energy,polyX45,polyY45,GNMPTS[14]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 16) { ingate = bantest(si[16].tpratio,si[16].energy,polyX46,polyY46,GNMPTS[15]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 17) { ingate = bantest(si[17].tpratio,si[17].energy,polyX47,polyY47,GNMPTS[16]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 18) { ingate = bantest(si[18].tpratio,si[18].energy,polyX48,polyY48,GNMPTS[17]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 19) { ingate = bantest(si[19].tpratio,si[19].energy,polyX49,polyY49,GNMPTS[18]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 20) { ingate = bantest(si[20].tpratio,si[20].energy,polyX410,polyY410,GNMPTS[19]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 21) { ingate = bantest(si[21].tpratio,si[21].energy,polyX411,polyY411,GNMPTS[20]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 22) { ingate = bantest(si[22].tpratio,si[22].energy,polyX412,polyY412,GNMPTS[21]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 23) { ingate = bantest(si[23].tpratio,si[23].energy,polyX413,polyY413,GNMPTS[22]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 24) { ingate = bantest(si[24].tpratio,si[24].energy,polyX414,polyY414,GNMPTS[23]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 25) { ingate = bantest(si[25].tpratio,si[25].energy,polyX415,polyY415,GNMPTS[24]);
                  if (ingate ==1)gaggvalid++;}
    if (i == 26) { ingate = bantest(si[26].tpratio,si[26].energy,polyX416,polyY416,GNMPTS[25]);
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

int gamproc(struct gdetector *ge, int ge_spe_xtl[MAX_GE_XTL*MAX_GE][8192], int ge_spe[MAX_GE][8192]){
/////////////////////
// G = Ge Detector //
/////////////////////      
int i,j,k;
long long int tdif;
int max1,max2,maxid1,maxid2;
int gmult,gcount;

      gmult = 0;
      for (i=0; i<MAX_GE; i++) {
        
        max1 = -1; max2 = -1;
        maxid1 = -1; maxid2 = -1;
        
        ge[i].id = i;
        ge[i].xvalid =0;  //C.W
      //Addback and Compton Suppression 
	      for (j=0; j<ge[i].xmult; j++) {
            
            //compton suppression per crystal 
             if (GE_BGO_SUPPRESSION == TRUE) {
                for (k=0; k<ge[i].bgomult; k++) {
                    tdif = abs( ge[i].xt[j] - ge[i].bgot[k] );
                    if (tdif >= 0 && tdif <= 4096) ge_bgo_tdif[i][tdif]++;
                    if ( (tdif < 50 && ge[i].bgoe[k] > 10 && ge[i].bgoe[k] < 30000) || ge[i].bgopileup[k]==TRUE) { //need to fix bgo pileup with trace analysis 
                        ge[i].xsuppress[j] = TRUE; 
                        ge[i].suppress = TRUE;                        
                    }  
                }
            }      
            
            //addback     
            if (ge[i].xsuppress[j] == FALSE) {       
                if (ge[i].xpileup[j] == TRUE) {
                    ge[i].pileup = TRUE;     
                    continue;
                }              
                //xtl spectra
                if (ge[i].xe[j] > 0 && ge[i].xe[j] < 8192 && ge[i].id >= 1) ge_spe_xtl[(ge[i].id-1)*MAX_GE_XTL + ge[i].xid[j]][ge[i].xe[j]]++;                        
                tdif = abs( ge[i].xt[j] - ge[i].time );
      //          if (tdif >= 0 && tdif <= 4096 && ge[i].time != 0 ) ge_xtl_tdif[i][tdif]++;           
                                       
                if (ge[i].xe[j] > 50 && ge[i].xe[j] < 5000) {                  
    //                if (tdif >= 0 && tdif <= 4096 && ge[i].time != 0 ) ge_xtl_tdif_ethresh[i][tdif]++;           
                    if (tdif < 20 || ge[i].time == 0) {  
                        ge[i].energy = ge[i].energy + ge[i].xe[j];
			                  ge[i].xevalid[ge[i].xvalid]=ge[i].xe[j]; //C.W 
                        //ge[i].thetavalid[xvalid]=ge[i].xtheta[j][1]; //C.W
                        ge[i].xidvalid[ge[i].xvalid]=ge[i].xid[j]; //C.W
                        ge[i].xvalid++; //C.W
                        ge[i].edop = ge[i].edop + (int) ge[i].xedopp[j]; //C.W
                        if (max1 < ge[i].xe[j]) {
                            max1 = ge[i].xe[j];
                            maxid1 = j;
                            ge[i].time = ge[i].xt[j];
                            ge[i].ctime = ge[i].xct[j];
                        } 
                    }
                    else {
                        ge[i].nonprompt = TRUE; // the first time will become the adopted value / event   
                    }        
                } 
            }//end for ge.xsuppress  --C. W  
           
	}//enf for germanium crystal multiplicity --C. W
        
        if (max1 == -1) continue;
        
        //Segmentation Position and Compton Suppression
        for (j=0; j<ge[i].smult; j++) {

            //compton suppression per segment 
            if (GE_BGO_SUPPRESSION == TRUE) {            
                for (k=0; k<ge[i].bgomult; k++) {
                    tdif = abs( ge[i].st[j] - ge[i].bgot[k] );
                    if (tdif < 50 && ge[i].bgoe[k] > 10) {
                        ge[i].ssuppress[j] = TRUE;    
                    }
                }
            }      

            //segment
            if (ge[i].ssuppress[j] == FALSE && ge[i].se[j] > 0 && ge[i].se[j] < 10000) {                    
                if (max2 < ge[i].se[j]) {
                    max2 = ge[i].se[j];
                    maxid2 = j;
                }                                           
            }  

        }    
        
        //Angle assignments
        ge[i].theta[0] = ge[i].xtheta[maxid1][0];                       //detector center
        ge[i].phi[0] = ge[i].xphi[maxid1][0];
        ge[i].theta[1] = ge[i].xtheta[maxid1][1];                       //crystal center
        ge[i].phi[1] = ge[i].xphi[maxid1][1];       
        if (ge[i].sid[maxid2] == 6) {                                   // side channel C  
            ge[i].theta[2] = ge[i].xtheta[maxid1][2];
            ge[i].phi[2] = ge[i].xphi[maxid1][2];  
        }
        else if (ge[i].sid[maxid2] == 5 || ge[i].sid[maxid2] == 7) {    // side channel L/R
            ge[i].theta[2] = ge[i].xtheta[maxid1][3];
            ge[i].phi[2] = ge[i].xphi[maxid1][3];             
        }  
        else {      
            ge[i].theta[2] = ge[i].xtheta[maxid1][1];                   // side channel failure --> crystal center
            ge[i].phi[2] = ge[i].xphi[maxid1][1];  
        }        
        
        //clean addback
//        if (ge[i].suppress == FALSE && ge[i].pileup == FALSE && ge[i].nonprompt == FALSE) {
        if (ge[i].suppress == FALSE && ge[i].nonprompt == FALSE) {

            ge[i].clean=TRUE; //maybe clean should not include suppress == FALSE?
        }    
        

        //Ge spectra
        if (ge[i].energy > 0 && ge[i].energy < 8192){ 
	    ge_spe[ge[i].id][ge[i].energy]++;
//	    ge_spe_dopp[ge[i].id][ge[i].edop]++; //C.W
	}         
  //      if (ge[i].energy > 0 && ge[i].energy < 4096 && ge[i].clean == TRUE) ge_spe_clean[ge[i].id][ge[i].energy]++;
        
        //copy data and increment counters
        ge[gmult]=ge[i];
        gmult++;        
        gcount++;
        
      } //end G
return gmult;
}



//Polarization Analysis:
void pol(int gaggvalid, int gmult, struct gdetector *ge, int pol_para[5000][5000], int pol_perp[5000][5000]){
if(gaggvalid == 1 && gmult >0){
int a,b;
for(b=0;b<gmult;b++){
	for(a=0;a<gmult;a++){
	if(b!=a){	 
       if(ge[a].xvalid == 2 && (ge[a].id == 3 || ge[a].id == 4  || ge[a].id == 9)){
                 //Parallel Spectra:
                  if((ge[a].xidvalid[0] == 1 && ge[a].xidvalid[1] == 4) || (ge[a].xidvalid[0] == 4 && ge[a].xidvalid[1] == 1)
                    || (ge[a].xidvalid[0] == 2 && ge[a].xidvalid[1] == 3) || (ge[a].xidvalid[0] == 3 && ge[a].xidvalid[1] == 2)){
                      ge[a].epol=ge[a].xevalid[0]+ge[a].xevalid[1];
                      printf("%f %d %d %d  %d\n",ge[a].theta[0],ge[a].xidvalid[0],ge[a].xidvalid[1],ge[a].xvalid,ge[a].epol);
                  //Histogram 1D spectra:
                   if((ge[a].epol > 0 && ge[a].epol < 5000) && (ge[b].energy > 0 && ge[b].energy < 5000))
                      pol_para[ge[b].energy][ge[a].epol]++;
                                                                   }
                
                  //Perpendicular Spectra:
                  if((ge[a].xidvalid[0] == 1 && ge[a].xidvalid[1] == 2) || (ge[a].xidvalid[0] == 2 && ge[a].xidvalid[1] == 1) 
                    || (ge[a].xidvalid[0] == 3 && ge[a].xidvalid[1] == 4) || (ge[a].xidvalid[0] == 4 && ge[a].xidvalid[1] == 3)){
                      ge[a].epol=ge[a].xevalid[0]+ge[a].xevalid[1];
                  //Histogram 1D spectra:
                      printf("%f %d %d %d  %d\n",ge[a].theta[0],ge[a].xidvalid[0],ge[a].xidvalid[1],ge[a].xvalid,ge[a].epol);
                   if((ge[a].epol > 0 && ge[a].epol <5000) && (ge[b].energy > 0 && ge[b].energy < 5000))
                      pol_perp[ge[b].energy][ge[a].epol]++;
                                                                  }
                                            }

              		}      
		}
	}

}
}


//Kinmeatic Correction:
void doppcorrect(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int option, float ml, float mbeam, float mh, float vbz, float betascale){
//Correcting Energy ://C.W
int gamrev;
if (gaggvalid == 1 && gmult > 0){
	//printf("%d %d\n",gaggvalid,gmult);
	for (gamrev=0;gamrev<gmult;gamrev++){
	//Start Doing Doppler Correction: //C.W Jun8:
	//Defining a Unit Vector:
	ge[gamrev].Dg[0]=sin((PIR)*ge[gamrev].theta[1])*cos((PIR)*ge[gamrev].phi[1]);
	ge[gamrev].Dg[1]=sin((PIR)*ge[gamrev].theta[1])*sin((PIR)*ge[gamrev].phi[1]);
	ge[gamrev].Dg[2]=cos((PIR)*ge[gamrev].theta[1]);

	//Deducing velocity of light particles:
	si[0].vl[0]=pow(2*(si[0].energy/100.)/ml,0.5)*sin((PIR)*si[0].theta[0])*cos((PIR)*si[0].phi[0]);
	si[0].vl[1]=pow(2*(si[0].energy/100.)/ml,0.5)*sin((PIR)*si[0].theta[0])*sin((PIR)*si[0].phi[0]);
	si[0].vl[2]=pow(2*(si[0].energy/100.)/ml,0.5)*cos((PIR)*si[0].theta[0]);


	//Deducing velocity of heavy particles:		
	si[0].vh[0]=-ml*si[0].vl[0]/mh;
	si[0].vh[1]=-ml*si[0].vl[1]/mh;
	si[0].vh[2]=(mbeam*vbz-ml*si[0].vl[2])/mh;

	si[0].betarest=pow(pow(si[0].vh[0],2)+pow(si[0].vh[1],2)+pow(si[0].vh[2],2),0.5);
	si[0].thetarest=100.0*acos((si[0].vh[2])/(si[0].betarest))*180.0/3.14159265;

	ge[gamrev].betavalue=(si[0].vh[0])*ge[gamrev].Dg[0]+(si[0].vh[1])*ge[gamrev].Dg[1]+(si[0].vh[2])*ge[gamrev].Dg[2];
	ge[gamrev].alpha=sin((PIR)*si[0].theta[0])*sin((PIR)*ge[gamrev].theta[1])*cos((PIR)*(si[0].phi[0]-ge[gamrev].phi[1]))+cos((PIR)*ge[gamrev].theta[1])*cos((PIR)*si[0].theta[0]);

//Non-relativistiv Doppler:
  if (option == 1) ge[gamrev].etrue=ge[gamrev].energy/(1+betascale*ge[gamrev].betavalue); //Doppler with kinmat
  else
	ge[gamrev].etrue=ge[gamrev].energy/(1+betascale*cos(PIR*ge[gamrev].theta[1])); //Doppler without Kinmat
//Relativistic Doppler Correction
	ge[gamrev].etrue2=ge[gamrev].energy*pow((1-pow(si[0].betarest,2)),0.5)/(1-(ge[gamrev].alpha)*(si[0].betarest));
//       printf("energync: %d,energyc: %f,betavalue1: %f, betavalue2: %f\n",ge[gamrev].energy,ge[gamrev].etrue,si[0].betarest,ge[gamrev].betavalue);

//Angle of Residual with respect to the beam axis:
						}			

printf("%s , energy: %d, etrue: %d\n","doppcorrect",ge[0].energy,(int)ge[0].etrue);

					}
}
//==========================================




void gammagagg1Dspectra(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int gagg_gamma_tdiff_ring2[4096][4096], int gagg_gamma_tdiff_ring4[4096][4096], int gagg_gamma_tdiff2[5000][5000],int gated_gamma[1][5000], int gated_gammatrue[1][5000], int gated_gamma_np[1][5000], int gated_gamma_nptrue[1][5000]){
int gem;
int gaggi;
//1 Dimensional Histogram //C. W
if (gaggvalid == 1  && gmult > 0){
	printf("energy: %d, etrue: %d\n",ge[0].energy, (int)ge[0].etrue);
        for (gem=0;gem<gmult;gem++){
           for (gaggi=0;gaggi<gaggvalid;gaggi++){
//put condition below if we want to have only one peak for the timediff spectra     
// if (si[gaggi].energy > 200 ){  
          si[gaggi].time=(si[gaggi].sit[0]+si[gaggi].sit[1])/2.;
         int tdiffgemgam=si[gaggi].time-ge[gem].time + 2000;
         if ((tdiffgemgam > 0 && tdiffgemgam < 4096)&&(si[gaggi].energy > 0 && si[gaggi].energy < 4096)){    
   	if (si[gaggi].id >= 1 && si[gaggi].id <=10)	
	    gagg_gamma_tdiff_ring2[(int)si[gaggi].energy][tdiffgemgam]++; 
	if (si[gaggi].id >=11 && si[gaggi].id <=26)
	    gagg_gamma_tdiff_ring4[(int)si[gaggi].energy][tdiffgemgam]++;
		}
         if ((tdiffgemgam > 0 && tdiffgemgam < 4096) && (ge[gem].energy > 0 && ge[gem].energy < 5000))
            gagg_gamma_tdiff2[ge[gem].energy][tdiffgemgam]++;
              
//Prompt gamma:
         // if ((ge[gem].energy > 0 && ge[gem].energy < 4096) && (tdiffgemgam > 1963 && tdiffgemgam < 2002)){
         if ((ge[gem].energy > 0 && ge[gem].energy < 5000) && (tdiffgemgam >=1960 && tdiffgemgam <=2002)){
	  gated_gamma[0][ge[gem].energy]++; 
          gated_gammatrue[0][(int)ge[gem].etrue]++;
	  ge[gem].validp = 1;
}
//Non-Prompt gamma:
       //   if ((ge[gem].energy > 0 && ge[gem].energy < 4096) &&( (tdiffgemgam > 2010 && tdiffgemgam < 2049))){
         if ((ge[gem].energy > 0 && ge[gem].energy < 5000) &&( (tdiffgemgam >= 1891 && tdiffgemgam <=1912) || (tdiffgemgam >= 2050 && tdiffgemgam <= 2071) )){
	  gated_gamma_np[0][ge[gem].energy]++;
	  gated_gamma_nptrue[0][(int)ge[gem].etrue]++;   
          ge[gem].validnp = 1;
            }                           
                            //  }
                    }
            }
}
}

void parthit(int gaggvalid, int gmult, struct gdetector *ge, struct sidetector *si, int e_lit2[MAX_SI][4096], int e_lit4[MAX_SI][4096], int theta_res2[MAX_SI][4096], int theta_res4[MAX_SI][4096], int gagghist[1][4096]){
//Histogram for Number of Hits:
if (gaggvalid >0 && gmult > 0){
if (gaggvalid == 1) {gagghist[0][1000]++;
if (si[0].id >= 1 && si[0].id <= 10){
	theta_res2[si[0].id][(int)si[0].thetarest]++;
	e_lit2[si[0].id][(int)si[0].energy]++;
//	printf("%d %d %d\n",gaggvalid,gmult,(int)si[0].thetarest);
}

//if (si[0].id >= 11) printf("%d\n",si[0].id);
if (si[0].id >= 11 && si[0].id <= 26){
	theta_res4[si[0].id][(int)si[0].thetarest]++;
	e_lit4[si[0].id][(int)si[0].energy]++;
//	printf("%d %d\n",gaggvalid,(int)si[0].thetarest);
}
}
if (gaggvalid == 2) gagghist[0][2000]++;
if (gaggvalid > 2) gagghist[0][3000]++;
}
}



