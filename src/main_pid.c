//Clarion2-Trinity Software Data Processing=====================================================
//C. Wibisono

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "global.h"
#include "param.h"
#include "spec.h"
#include "spec_dec.h"
#include "gaggpid.h"
#include "gaggpid_dec.h"

struct subevent subevt[MAX_ID]={0};
int sevtmult=0;
unsigned long long int sevtcount=0;
unsigned long long int pileupcount=0;
unsigned long long int evtcount=0;


struct gdetector ge[MAX_GE]={0};
int gmult=0;
unsigned long long int gcount=0;



struct sidetector si[MAX_SI]={0};
//int simult=0;
unsigned long long int sicount=0;
//====================================================

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
  float Ebeam,mb,mbeam,ml0,ml,mh0,mh,amu,vbz;
  int parttype=0;
  amu=931.5;
 /*
  Ebeam=30;
  mbeam=16*amu;
  ml=1*amu;
  mh=33*amu;
*/
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
  //Allocate Memory for PID making:
  int ab;
  for(ab=0;ab<4096;ab++){
	if((pid_qdc21[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc22[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc23[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc24[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc25[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc26[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc27[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc28[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc29[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc210[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc41[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc42[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc43[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc44[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc45[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc46[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc47[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc48[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc49[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc410[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc411[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc412[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc413[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc414[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc415[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
	if((pid_qdc416[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}
   }//end iteration for malloc
 
   
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
FILE *fgagg, *fparticle;
fgagg=fopen(argv[5],"r");
int idgagg;
int n;
float quadgagg,slopegagg,interceptgagg;
//int parttype; //define particle type
char type[20];

fparticle=fopen(argv[9],"r");

if (fparticle == NULL){
fprintf(stderr,"Error, cannot open input file %s\n",argv[9]);
}

while(fgets(line,LINE_LENGTH,fparticle) !=NULL){
        if(strchr(line,'#') ==NULL)
                sscanf(line,"%f\t%f\t%f\t%f\t%s\t%d\n",&Ebeam,&mb,&ml0,&mh0,type,&parttype);
                        }
memset(line,0,LINE_LENGTH);
mbeam=mb*amu; //conversion to amu
ml=ml0*amu;
mh=mh0*amu;
vbz=pow(2*Ebeam/mbeam,0.5);

if (fgagg == NULL){
fprintf(stderr,"Error, cannot open input file %s\n",argv[5]);
}

//Protons:
if (parttype == 1){
while (fgets(line,LINE_LENGTH,fgagg) != NULL){
sscanf(line,"%d\t%f\t%f\n",&idgagg,&slopegagg,&interceptgagg);
gaid[n]=idgagg;
gaggslope[n]=slopegagg;
gaggintercept[n]=interceptgagg;
n++;
memset(line,0,LINE_LENGTH);
	}
}


//Alphas:
if (parttype == 2){
while (fgets(line,LINE_LENGTH,fgagg) != NULL){
sscanf(line,"%d\t%f\t%f\t%f\n",&idgagg,&interceptgagg,&slopegagg,&quadgagg);
gaid[n]=idgagg;
gaggquad[n]=quadgagg;
gaggslope[n]=slopegagg;
gaggintercept[n]=interceptgagg;
n++;
memset(line,0,LINE_LENGTH);
}
}


fclose(fgagg);
fclose(fparticle);
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


//Read to file and build event:
sevtmult=toreader(sub,subevt,fpr);

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
//printf("%d\n",parttype);

//Create the PID:
gaggpid(si,parttype);
/*
//Checking Gagg and Germanium Multiplicities:
int gaggvalid=gaggproc(si,parttype);
gmult=gamproc(ge);



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

*/

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

//Perform Polarization Analysis:      
//pol(gaggvalid, gmult,ge,pol_para, pol_perp);

//Perform DCO Ratio:
//dco(gmult,gaggvalid,ge);

//Perform DCO Kinmat:
//dcokinmat(gmult,gaggvalid,ge);
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
/*
     write_data4(GE_SPE_CLEAN, *ge_spe_clean, MAX_GE, 4096, overwrite);  
    write_data4(PID_Gamma_R, *ge_gagg_spe_r, MAX_GE, 4096, overwrite);
    write_data4(PID_Gammaa, *gagg_gamma_tdiff_ring2, 4096, 4096, overwrite);
    write_data4(PID_Gammab, *gagg_gamma_tdiff_ring4, 4096, 4096, overwrite);
    write_data4(PID_Gamma2, *gagg_gamma_tdiff2, 5000, 5000, overwrite);
    write_data4(PID_Gamma_1Dp, *gated_gamma, 1, 5000, overwrite);
    write_data4(PID_Gamma_1Dptrue, *gated_gammatrue,1,5000,overwrite);
    write_data4(PID_Gamma_1Dnp, *gated_gamma_np, 1, 5000, overwrite);
    write_data4(PID_Gamma_1Dnptrue, *gated_gamma_nptrue, 1, 5000, overwrite);
*/
/*
  for (int i=0; i<26; ++i) {
    char name[128];
    sprintf(name, "gagg_pid_%i.spn", i+1);
    write_data4(name, &pid_all[i][0][0], 4096, 4096, overwrite);
  }

*/
 
 
 write_data4dyn(PID_QDC21, pid_qdc21, 4096,4096, overwrite);
 write_data4dyn(PID_QDC22, pid_qdc22, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC23, pid_qdc23, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC24,pid_qdc24, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC25, pid_qdc25, 4096,4096, overwrite); 
 write_data4dyn(PID_QDC26, pid_qdc26, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC27, pid_qdc27, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC28,pid_qdc28, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC29, pid_qdc29, 4096,4096, overwrite);
 write_data4dyn(PID_QDC210, pid_qdc210, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC41, pid_qdc41, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC42,pid_qdc42, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC43, pid_qdc43, 4096,4096, overwrite);
 write_data4dyn(PID_QDC44, pid_qdc44, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC45, pid_qdc45, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC46,pid_qdc46, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC47, pid_qdc47, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC48, pid_qdc48, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC49, pid_qdc49, 4096,4096, overwrite);
 write_data4dyn(PID_QDC410, pid_qdc410, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC411, pid_qdc411, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC412, pid_qdc412, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC413, pid_qdc413, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC414, pid_qdc414, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC415, pid_qdc415, 4096, 4096, overwrite);
 write_data4dyn(PID_QDC416, pid_qdc416, 4096, 4096, overwrite);

// write_data4(GAGG_GE_TDIFF, *gagg_ge_tdiff, 4096, 4096, overwrite);

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
 /*
  write_data4(GAGG_HIST, *gagghist, 4096, 1, overwrite);
  write_data4(THETA_RES2, *theta_res2, 4096, 28, overwrite);
  write_data4(THETA_RES4, *theta_res4, 4096, 28, overwrite);
  write_data4(E_LIT2, *e_lit2, 4096, 28, overwrite);
  write_data4(E_LIT4, *e_lit4, 4096, 28, overwrite);
  write_data4(POL_PARA,*pol_para, 5000, 5000, overwrite);
  write_data4(POL_PERP, *pol_perp, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PUREF, *gg_prompt_puref, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PUREFP, *gg_prompt_purefp, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PUREFPKIN, *gg_prompt_purefpkin, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PUREB, *gg_prompt_pureb, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PUREBP, *gg_prompt_purebp, 5000, 5000, overwrite);
  write_data4(GG_PROMPT_PUREBPKIN, *gg_prompt_purebpkin, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PUREF, *gg_nprompt_puref, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PUREFP, *gg_nprompt_purefp, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PUREFPKIN, *gg_nprompt_purefpkin, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PUREB, *gg_nprompt_pureb, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PUREBP, *gg_nprompt_purebp, 5000, 5000, overwrite);
  write_data4(GG_NPROMPT_PUREBPKIN, *gg_nprompt_purebpkin, 5000, 5000, overwrite);
*/
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






