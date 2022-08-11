//Write Histogram Function:
//Refurbished from ORNL Program

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




