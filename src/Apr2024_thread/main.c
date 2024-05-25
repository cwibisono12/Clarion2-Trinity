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
#include "ggmat.h"
#include "ggmatglobal.h"

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
  int ab,bc;
  
  for(ab=0;ab<4096;ab++){
	if((gg_tdif1[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_tdif2[ab]=(int*)malloc(4096*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	}//end iteration for malloc
 
  for(bc=0;bc<5000;bc++){
	if((gg_prompt[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_prompt_nddopp[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_prompt_ddopp[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt1[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_nddopp1[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_ddopp1[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt2[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_nddopp2[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_ddopp2[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_nddopp[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_ddopp[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_prompt_ng[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_prompt_nddopp_ng[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_prompt_ddopp_ng[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_ng[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_nddopp_ng[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
	printf("\nError,memory not allocated\n");
	exit(1);
	}   
	if((gg_nprompt_ddopp_ng[bc]=(int*)malloc(5000*sizeof(int))) == NULL){
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
sscanf(line,"%d\t%f\t%f\t%f\n",&idgagg,&quadgagg,&slopegagg,&interceptgagg);
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
    g22=0,g23=0,g24=0,g25=0,g26=0,g27=0,g28=0,g29=0,g30=0,
    g31=0,g32=0,g33=0,g34=0,g35=0,g36=0,g37=0,g38=0,g39=0,g40=0,
    g41=0,g42=0,g43=0,g44=0,g45=0,g46=0,g47=0,g48=0,g49=0,g50=0,
    g51=0,g52=0,g53=0,g54=0,g55=0,g56=0,g57=0,g58=0,g59=0,g60=0,
    g61=0,g62=0,g63=0,g64=0;

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
        polyX11[g1]=banx;
        polyY11[g1]=bany;
      printf("%f %f\n",polyX11[g1],polyY11[g1]); 
       g1++;
        }
    if (gaggid == GID[1] && gaggnmpts == GNMPTS[1]){
        polyX12[g2]=banx;
        polyY12[g2]=bany;
      printf("%f %f\n",polyX12[g2],polyY12[g2]); 
       g2++;
         }
    if (gaggid == GID[2] && gaggnmpts == GNMPTS[2]){
        polyX13[g3]=banx;
        polyY13[g3]=bany;
       printf("%f %f\n",polyX13[g3],polyY13[g3]);
        g3++;
          }
    if (gaggid == GID[3] && gaggnmpts == GNMPTS[3]){
        polyX14[g4]=banx;
        polyY14[g4]=bany;
        printf("%f %f\n",polyX14[g4],polyY14[g4]);
        g4++;
          }
    if (gaggid == GID[4] && gaggnmpts == GNMPTS[4]){
        polyX15[g5]=banx;
        polyY15[g5]=bany;
        printf("%f %f\n",polyX15[g5],polyY15[g5]);
        g5++;
          }
    if (gaggid == GID[5] && gaggnmpts == GNMPTS[5]){
        polyX16[g6]=banx;
        polyY16[g6]=bany;
        printf("%f %f\n",polyX16[g6],polyY16[g6]);
        g6++;
          }
    if (gaggid == GID[6] && gaggnmpts == GNMPTS[6]){
        polyX17[g7]=banx;
        polyY17[g7]=bany;
        printf("%f %f\n",polyX17[g7],polyY17[g7]);
        g7++;
          }
    if (gaggid == GID[7] && gaggnmpts == GNMPTS[7]){
        polyX18[g8]=banx;
        polyY18[g8]=bany;
        printf("%f %f\n",polyX18[g8],polyY18[g8]);
        g8++;
          }
    if (gaggid == GID[8] && gaggnmpts == GNMPTS[8]){
        polyX21[g9]=banx;
        polyY21[g9]=bany;
        printf("%f %f\n",polyX21[g9],polyY21[g9]); 
       g9++;
        }
    if (gaggid == GID[9] && gaggnmpts == GNMPTS[9]){
        polyX22[g10]=banx;
        polyY22[g10]=bany;
      printf("%f %f\n",polyX22[g10],polyY22[g10]); 
       g10++;
         }
    if (gaggid == GID[10] && gaggnmpts == GNMPTS[10]){
        polyX23[g11]=banx;
        polyY23[g11]=bany;
       printf("%f %f\n",polyX23[g11],polyY23[g11]);
        g11++;
          }
    if (gaggid == GID[11] && gaggnmpts == GNMPTS[11]){
        polyX24[g12]=banx;
        polyY24[g12]=bany;
        printf("%f %f\n",polyX24[g12],polyY14[g12]);
        g12++;
          }
    if (gaggid == GID[12] && gaggnmpts == GNMPTS[12]){
        polyX25[g13]=banx;
        polyY25[g13]=bany;
        printf("%f %f\n",polyX25[g13],polyY25[g13]);
        g13++;
          }
    if (gaggid == GID[13] && gaggnmpts == GNMPTS[13]){
        polyX26[g14]=banx;
        polyY26[g14]=bany;
        printf("%f %f\n",polyX26[g14],polyY26[g14]);
        g14++;
          }
    if (gaggid == GID[14] && gaggnmpts == GNMPTS[14]){
        polyX27[g15]=banx;
        polyY27[g15]=bany;
        printf("%f %f\n",polyX27[g15],polyY27[g15]);
        g15++;
          }
    if (gaggid == GID[15] && gaggnmpts == GNMPTS[15]){
        polyX28[g16]=banx;
        polyY28[g16]=bany;
        printf("%f %f\n",polyX28[g16],polyY28[g16]);
        g16++;
          }
    if (gaggid == GID[16] && gaggnmpts == GNMPTS[16]){
        polyX29[g17]=banx;
        polyY29[g17]=bany;
        printf("%f %f\n",polyX29[g17],polyY29[g17]);
        g17++;
          }
    if (gaggid == GID[17] && gaggnmpts == GNMPTS[17]){
        polyX210[g18]=banx;
        polyY210[g18]=bany;
        printf("%f %f\n",polyX210[g18],polyY210[g18]);
        g18++;
          }
    if (gaggid == GID[18] && gaggnmpts == GNMPTS[18]){
        polyX31[g19]=banx;
        polyY31[g19]=bany;
        printf("%f %f\n",polyX31[g19],polyY210[g19]);
        g19++;
          }
    if (gaggid == GID[19] && gaggnmpts == GNMPTS[19]){
        polyX32[g20]=banx;
        polyY32[g20]=bany;
        printf("%f %f\n",polyX32[g20],polyY32[g20]);
        g20++;
          }
    if (gaggid == GID[20] && gaggnmpts == GNMPTS[20]){
        polyX33[g21]=banx;
        polyY33[g21]=bany;
        printf("%f %f\n",polyX33[g21],polyY33[g21]);
        g21++;
          }
    if (gaggid == GID[21] && gaggnmpts == GNMPTS[21]){
        polyX34[g22]=banx;
        polyY34[g22]=bany;
        printf("%f %f\n",polyX34[g22],polyY34[g22]);
        g22++;
          }
    if (gaggid == GID[22] && gaggnmpts == GNMPTS[22]){
        polyX35[g23]=banx;
        polyY35[g23]=bany;
        printf("%f %f\n",polyX35[g23],polyY35[g23]);
        g23++;
          }
    if (gaggid == GID[23] && gaggnmpts == GNMPTS[23]){
        polyX36[g24]=banx;
        polyY36[g24]=bany;
        printf("%f %f\n",polyX36[g24],polyY32[g24]);
        g24++;
          }
    if (gaggid == GID[24] && gaggnmpts == GNMPTS[24]){
        polyX37[g25]=banx;
        polyY37[g25]=bany;
        printf("%f %f\n",polyX37[g25],polyY37[g25]);
        g25++;
          }
    if (gaggid == GID[25] && gaggnmpts == GNMPTS[25]){
        polyX38[g26]=banx;
        polyY38[g26]=bany;
        printf("%f %f\n",polyX38[g26],polyY38[g26]);
        g26++;
          }
    if (gaggid == GID[26] && gaggnmpts == GNMPTS[26]){
        polyX39[g27]=banx;
        polyY39[g27]=bany;
        printf("%f %f\n",polyX39[g27],polyY39[g27]);
        g27++;
          }
    if (gaggid == GID[27] && gaggnmpts == GNMPTS[27]){
        polyX310[g28]=banx;
        polyY310[g28]=bany;
        printf("%f %f\n",polyX310[g28],polyY310[g28]);
        g28++;
          }
    if (gaggid == GID[28] && gaggnmpts == GNMPTS[28]){
        polyX311[g29]=banx;
        polyY311[g29]=bany;
        printf("%f %f\n",polyX311[g29],polyY311[g29]);
        g29++;
          }
    if (gaggid == GID[29] && gaggnmpts == GNMPTS[29]){
        polyX312[g30]=banx;
        polyY312[g30]=bany;
        printf("%f %f\n",polyX312[g30],polyY312[g30]);
        g30++;
          }
    if (gaggid == GID[30] && gaggnmpts == GNMPTS[30]){
        polyX313[g31]=banx;
        polyY313[g31]=bany;
        printf("%f %f\n",polyX313[g31],polyY313[g31]);
        g31++;
          }
    if (gaggid == GID[31] && gaggnmpts == GNMPTS[31]){
        polyX314[g32]=banx;
        polyY314[g32]=bany;
        printf("%f %f\n",polyX314[g32],polyY314[g32]);
        g32++;
          }
    if (gaggid == GID[32] && gaggnmpts == GNMPTS[32]){
        polyX41[g33]=banx;
        polyY41[g33]=bany;
        printf("%f %f\n",polyX41[g33],polyY41[g33]);
        g33++;
          }
    if (gaggid == GID[33] && gaggnmpts == GNMPTS[33]){
        polyX42[g34]=banx;
        polyY42[g34]=bany;
        printf("%f %f\n",polyX42[g34],polyY42[g34]);
        g34++;
          }
    if (gaggid == GID[34] && gaggnmpts == GNMPTS[34]){
        polyX43[g35]=banx;
        polyY43[g35]=bany;
        printf("%f %f\n",polyX43[g35],polyY43[g35]);
        g35++;
          }
    if (gaggid == GID[35] && gaggnmpts == GNMPTS[35]){
        polyX44[g36]=banx;
        polyY44[g36]=bany;
        printf("%f %f\n",polyX44[g36],polyY44[g36]);
        g36++;
          }
    if (gaggid == GID[36] && gaggnmpts == GNMPTS[36]){
        polyX45[g37]=banx;
        polyY45[g37]=bany;
        printf("%f %f\n",polyX45[g37],polyY45[g37]);
        g37++;
          }
    if (gaggid == GID[37] && gaggnmpts == GNMPTS[37]){
        polyX46[g38]=banx;
        polyY46[g38]=bany;
        printf("%f %f\n",polyX46[g38],polyY46[g38]);
        g38++;
          }
    if (gaggid == GID[38] && gaggnmpts == GNMPTS[38]){
        polyX47[g39]=banx;
        polyY47[g39]=bany;
        printf("%f %f\n",polyX47[g39],polyY47[g39]);
        g39++;
          }
    if (gaggid == GID[39] && gaggnmpts == GNMPTS[39]){
        polyX48[g40]=banx;
        polyY48[g40]=bany;
        printf("%f %f\n",polyX48[g40],polyY48[g40]);
        g40++;
          }
    if (gaggid == GID[40] && gaggnmpts == GNMPTS[40]){
        polyX49[g41]=banx;
        polyY49[g41]=bany;
        printf("%f %f\n",polyX49[g41],polyY49[g41]);
        g41++;
          }
    if (gaggid == GID[41] && gaggnmpts == GNMPTS[41]){
        polyX410[g42]=banx;
        polyY410[g42]=bany;
        printf("%f %f\n",polyX410[g42],polyY410[g42]);
        g42++;
          }
    if (gaggid == GID[42] && gaggnmpts == GNMPTS[42]){
        polyX411[g43]=banx;
        polyY411[g43]=bany;
        printf("%f %f\n",polyX411[g43],polyY411[g43]);
        g43++;
          }
    if (gaggid == GID[43] && gaggnmpts == GNMPTS[43]){
        polyX412[g44]=banx;
        polyY412[g44]=bany;
        printf("%f %f\n",polyX412[g44],polyY412[g44]);
        g44++;
          }          
    if (gaggid == GID[44] && gaggnmpts == GNMPTS[44]){
        polyX413[g45]=banx;
        polyY413[g45]=bany;
        printf("%f %f\n",polyX413[g45],polyY413[g45]);
        g45++;
          }
    if (gaggid == GID[45] && gaggnmpts == GNMPTS[45]){
        polyX414[g46]=banx;
        polyY414[g46]=bany;
        printf("%f %f\n",polyX414[g46],polyY414[g46]);
        g46++;
          }
    if (gaggid == GID[46] && gaggnmpts == GNMPTS[46]){
        polyX415[g47]=banx;
        polyY415[g47]=bany;
        printf("%f %f\n",polyX415[g47],polyY415[g47]);
        g47++;
          }
    if (gaggid == GID[47] && gaggnmpts == GNMPTS[47]){
        polyX416[g48]=banx;
        polyY416[g48]=bany;
        printf("%f %f\n",polyX416[g48],polyY416[g48]);
        g48++;
          }
    if (gaggid == GID[48] && gaggnmpts == GNMPTS[48]){
        polyX51[g49]=banx;
        polyY51[g49]=bany;
        printf("%f %f\n",polyX51[g49],polyY51[g49]);
        g49++;
          }
    if (gaggid == GID[49] && gaggnmpts == GNMPTS[49]){
        polyX52[g50]=banx;
        polyY52[g50]=bany;
        printf("%f %f\n",polyX52[g50],polyY52[g50]);
        g50++;
          }
    if (gaggid == GID[50] && gaggnmpts == GNMPTS[50]){
        polyX53[g51]=banx;
        polyY53[g51]=bany;
        printf("%f %f\n",polyX53[g51],polyY53[g51]);
        g51++;
          }
    if (gaggid == GID[51] && gaggnmpts == GNMPTS[51]){
        polyX54[g52]=banx;
        polyY54[g52]=bany;
        printf("%f %f\n",polyX54[g52],polyY54[g52]);
        g52++;
          }
    if (gaggid == GID[52] && gaggnmpts == GNMPTS[52]){
        polyX55[g53]=banx;
        polyY55[g53]=bany;
        printf("%f %f\n",polyX55[g53],polyY55[g53]);
        g53++;
          }
    if (gaggid == GID[53] && gaggnmpts == GNMPTS[53]){
        polyX56[g54]=banx;
        polyY56[g54]=bany;
        printf("%f %f\n",polyX56[g54],polyY56[g54]);
        g54++;
          }
    if (gaggid == GID[54] && gaggnmpts == GNMPTS[54]){
        polyX57[g55]=banx;
        polyY57[g55]=bany;
        printf("%f %f\n",polyX57[g55],polyY57[g55]);
        g55++;
          }
    if (gaggid == GID[55] && gaggnmpts == GNMPTS[55]){
        polyX58[g56]=banx;
        polyY58[g56]=bany;
        printf("%f %f\n",polyX58[g56],polyY58[g56]);
        g56++;
          }
    if (gaggid == GID[56] && gaggnmpts == GNMPTS[56]){
        polyX59[g57]=banx;
        polyY59[g57]=bany;
        printf("%f %f\n",polyX59[g57],polyY59[g57]);
        g57++;
          }
    if (gaggid == GID[57] && gaggnmpts == GNMPTS[57]){
        polyX510[g58]=banx;
        polyY510[g58]=bany;
        printf("%f %f\n",polyX510[g58],polyY510[g58]);
        g58++;
          }
    if (gaggid == GID[58] && gaggnmpts == GNMPTS[58]){
        polyX511[g59]=banx;
        polyY511[g59]=bany;
        printf("%f %f\n",polyX511[g59],polyY511[g59]);
        g59++;
          }
    if (gaggid == GID[59] && gaggnmpts == GNMPTS[59]){
        polyX512[g60]=banx;
        polyY512[g60]=bany;
        printf("%f %f\n",polyX512[g60],polyY512[g60]);
        g60++;
          }          
    if (gaggid == GID[60] && gaggnmpts == GNMPTS[60]){
        polyX513[g61]=banx;
        polyY513[g61]=bany;
        printf("%f %f\n",polyX513[g61],polyY513[g61]);
        g61++;
          }
    if (gaggid == GID[61] && gaggnmpts == GNMPTS[61]){
        polyX514[g62]=banx;
        polyY514[g62]=bany;
        printf("%f %f\n",polyX514[g62],polyY514[g62]);
        g62++;
          }
    if (gaggid == GID[62] && gaggnmpts == GNMPTS[62]){
        polyX515[g63]=banx;
        polyY515[g63]=bany;
        printf("%f %f\n",polyX515[g63],polyY515[g63]);
        g63++;
          }
    if (gaggid == GID[63] && gaggnmpts == GNMPTS[63]){
        polyX516[g64]=banx;
        polyY516[g64]=bany;
        printf("%f %f\n",polyX516[g64],polyY516[g64]);
        g64++;
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
//gaggpid(si,parttype);

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
//parthit(gaggvalid, gmult, ge, si);


//Perform 1D gagg-Gamma Histogram:
gammagagg1Dspectra(gaggvalid, gmult, si, ge);

//Make gg matrix:
ggmat(gaggvalid,gmult,ge);

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
  */
  write_data4(E_RAW, *e_raw, MAX_ID, 8192, overwrite);
  write_data4(E_CAL, *e_cal, MAX_ID, 8192, overwrite);
/*
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
 
  write_data4(GE_BGO_TDIF, *ge_bgo_tdif, MAX_GE, 4096, overwrite);
  write_data4(GE_XTL_TDIF, *ge_xtl_tdif, MAX_GE, 4096, overwrite);
  write_data4(GE_XTL_TDIF_ETHRESH, *ge_xtl_tdif_ethresh, MAX_GE, 4096, overwrite);
  write_data4(GE_SPE_XTL, *ge_spe_xtl, MAX_GE*MAX_GE_XTL, 8192, overwrite);
  write_data4(GE_SPE, *ge_spe, MAX_GE, 8192, overwrite);
  write_data4(GE_SPE_DOPP, *ge_spe_dopp, MAX_GE, 8192, overwrite); //C.W
  write_data4(GE_SPE_CLEAN, *ge_spe_clean, MAX_GE, 8192, overwrite);

/* 
 //trinity
  write_data4(PID, *pid, 4096, 4096, overwrite);
  write_data4(PID_EVSP, *pid_evsp, 4096, 4096, overwrite);
  write_data4(PID_EVST, *pid_evst, 4096, 4096, overwrite);
  write_data4(PID_EVSTAU, *pid_evstau, 4096, 4096, overwrite);
  write_data4(PID_EVSR, *pid_evsr, 4096, 4096, overwrite);
  write_data4(PID_TAUVSR, *pid_tauvsr, 4096, 4096, overwrite);
  write_data4(GAGGDT, *gaggdt, 4096, 4096, overwrite);
*/

    //write_data4(GE_SPE_CLEAN, *ge_spe_clean, MAX_GE, 4096, overwrite);  
   // write_data4(PID_Gamma_R, *ge_gagg_spe_r, MAX_GE, 4096, overwrite);
    //write_data4(GAGG_GAMMA_TDIFF_RING2, *gagg_gamma_tdiff_ring2, 4096, 4096, overwrite);
    //write_data4(GAGG_GAMMA_TDIFF_RING4, *gagg_gamma_tdiff_ring4, 4096, 4096, overwrite);
    write_data4(GAGG_GAMMA_TDIFF2, *gagg_gamma_tdiff2, 5000, 5000, overwrite);
    write_data4(GATED_GAMMA, *gated_gamma, 16, 5000, overwrite);
    write_data4(GATED_GAMMATRUE, *gated_gammatrue,16,5000,overwrite);
    write_data4(GATED_GAMMA_NP, *gated_gamma_np, 16, 5000, overwrite);
    write_data4(GATED_GAMMA_NPTRUE, *gated_gamma_nptrue, 16, 5000, overwrite);

    write_data4(GATED_GAMMARING,*gated_gammaRing,66,5000,overwrite);
    write_data4(GATED_GAMMARINGTRUE,*gated_gammaRingtrue,66,5000,overwrite);
    write_data4(GATED_GAMMARING_NP,*gated_gammaRing_np,66,5000,overwrite);
    write_data4(GATED_GAMMARING_NPTRUE,*gated_gammaRing_nptrue,66,5000,overwrite);

/*
  for (int i=0; i<26; ++i) {
    char name[128];
    sprintf(name, "gagg_pid_%i.spn", i+1);
    write_data4(name, &pid_all[i][0][0], 4096, 4096, overwrite);
  }

*/
 
 

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

 // write_data4dyn(GG_TDIF, gg_tdif, 4096, 4096, overwrite);
  write_data4dyn(GG_TDIF1, gg_tdif1, 4096, 4096, overwrite);
  write_data4dyn(GG_TDIF2, gg_tdif2, 4096, 4096, overwrite);
  write_data4dyn(GG_PROMPT, gg_prompt, 5000, 5000, overwrite);
 // write_data4dyn(GG_PROMPT_PURE, gg_prompt_pure, 5000, 5000, overwrite);
 // write_data4dyn(GG_NPROMPT_PURE, gg_nprompt_pure, 5000, 5000, overwrite);
  write_data4dyn(GG_NPROMPT, gg_nprompt, 5000, 5000, overwrite);
  write_data4dyn(GG_NPROMPT1, gg_nprompt1, 5000, 5000, overwrite);
  write_data4dyn(GG_NPROMPT2, gg_nprompt2, 5000, 5000, overwrite);
  write_data4dyn(GG_PROMPT_NDDOPP, gg_prompt_nddopp, 5000, 5000, overwrite); //C.W
  write_data4dyn(GG_PROMPT_DDOPP, gg_prompt_ddopp, 5000, 5000, overwrite); //C. W
  write_data4dyn(GG_NPROMPT_NDDOPP, gg_nprompt_nddopp, 5000,5000, overwrite); //C.W
  write_data4dyn(GG_NPROMPT_NDDOPP1, gg_nprompt_nddopp1, 5000,5000, overwrite); //C.W
  write_data4dyn(GG_NPROMPT_DDOPP1, gg_nprompt_ddopp1, 5000, 5000, overwrite); //C.W
  write_data4dyn(GG_NPROMPT_NDDOPP2, gg_nprompt_nddopp2, 5000,5000, overwrite); //C.W
  write_data4dyn(GG_NPROMPT_DDOPP2,gg_nprompt_ddopp2, 5000, 5000, overwrite); //C.W
  write_data4dyn(GG_NPROMPT_DDOPP,gg_nprompt_ddopp, 5000, 5000, overwrite);//C.W
  write_data4dyn(GG_PROMPT_NG, gg_prompt_ng, 5000, 5000, overwrite);
  write_data4dyn(GG_PROMPT_NDDOPP_NG, gg_prompt_nddopp_ng, 5000, 5000, overwrite); //C.W
  write_data4dyn(GG_PROMPT_DDOPP_NG, gg_prompt_ddopp_ng, 5000, 5000, overwrite); //C. W
  write_data4dyn(GG_NPROMPT_NG, gg_nprompt_ng, 5000, 5000, overwrite);
  write_data4dyn(GG_NPROMPT_NDDOPP_NG, gg_nprompt_nddopp_ng, 5000,5000, overwrite); //C.W
  write_data4dyn(GG_NPROMPT_DDOPP_NG,gg_nprompt_ddopp_ng, 5000, 5000, overwrite);//C.W
 
  fclose(fpr);
  fclose(debugfile);
  
  return 0;
}






