//Write Histogram Function:
//write_data2 and write_data4 are based on work J. M Allmond
//write_data4dyn is based on work C. Wibisono


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void write_data4dyn(char *filename, int *data[], int xdim, int ydim, int overwrite){
    FILE *FP;
    int i,j;
    int *previous[ydim];
   
    if(!overwrite) {
	for(i=0;i<ydim;i++){
        //allocate memory for 1d-array for reading in rows of 2d Radware matrix
        if ( ( previous[i] = (int *)malloc(xdim *  sizeof(int)) ) == NULL ) {
            printf("\nError, memory not allocated.\n");
            exit(1);
        }
      	}
        //open previous spectra file  
        if( (FP=fopen(filename, "r")) != NULL ){
           for(i=0;i<ydim;i++){ 
	   fread(previous[i], sizeof(int)*xdim, 1, FP);        
            }
	    fclose(FP);
            //update spectra
       for(j=0;j<xdim;j++){    
	for (i=0; i<ydim; i++) {
                if(previous[i][j] < (powf(2,sizeof(int)*8.0)-2)) 
                    data[i][j] = data[i][j] + previous[i][j];
            }   
        } 
	}
        else{
            printf("%s did not previously exist, creating ...\n", filename);       
        }   
	//printf("hello memory\n"); 
      
        //Deallocate previous data
      for(i=0;i<ydim;i++){
//	 if(previous[i] != NULL) 
	 free(previous[i]);
	}
    }
  //printf("hello memory2\n");
    FP=fopen(filename, "w");
    for(i=0;i<ydim;i++){
    fwrite(data[i], sizeof(int)*xdim,1 , FP);
    }
    fclose(FP);    
}



