//List of Declaration for Gloval Variables//
//C. Wibisono
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
void write_data4dyn(char *filename, int *data[], int xdim, int ydim, int overwrite); //C.W

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


int toreader(unsigned int sub[], struct subevent *subevt, FILE *fpr);

void detmaps(int sevtmult, struct subevent *subevt, struct gdetector *ge, struct sidetector *si);

int gaggproc(struct sidetector *si,int parttype);

int gamproc(struct gdetector *ge);

void doppcorrect(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge, int option, float ml, float mbeam, float mh, float vbz, float betascale);

void gammagagg1Dspectra(int gaggvalid, int gmult, struct sidetector *si, struct gdetector *ge);

void parthit(int gaggvalid, int gmult, struct gdetector *ge, struct sidetector *si);

void ggmat(int gaggvalid, int gmult, struct gdetector *ge);







