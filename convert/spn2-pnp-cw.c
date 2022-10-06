/***********************************************************/
/*  spn-pnp: J.M. Allmond //revised_Version: cw19s         */
/*                                                         */
/*  subtract non-prompt from prompt with a scale factor    */
/*                                                         */
/*  Output -->    output.spn                               */
/*                                                         */
/*  Compile -->   gcc -o spn-pnp -lm spn-pnp.c             */
/*
/RUn: ./spn-pnp-cw prompt nprompt nprompt1 nprompt2 1.0                                                         */ 
/***********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#define MAT_TYPE int 

#define Y_DIM 5000 
#define X_DIM 5000


MAT_TYPE array2d1[Y_DIM][X_DIM]={{0}};
MAT_TYPE array2d2[Y_DIM][X_DIM]={{0}};
MAT_TYPE array2d3[Y_DIM][X_DIM]={{0}};
MAT_TYPE array2d4[Y_DIM][X_DIM]={{0}};
MAT_TYPE array2d5[Y_DIM][X_DIM]={{0}};


/////////////////////////////////////////////////////////////
// START MAIN PROGRAM                                      //    
/////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
    
FILE *fp;
                        
int i,j,k;
float scale;

if (argc==6) {
    scale = atof(argv[5]);
}
else {
    printf("./spn-pnp prompt.spn nonprompt.spn scalefactor\n");
    exit(1);    
}    

if( (fp=fopen(argv[1], "r")) == NULL ){
    printf("\nError, could not open input file.\n");    
    exit(1);    
}    
fread(array2d1, sizeof(array2d1), 1, fp); 
fclose(fp);

if( (fp=fopen(argv[2], "r")) == NULL ){
    printf("\nError, could not open input file.\n");    
    exit(1);    
}    
fread(array2d2, sizeof(array2d2), 1, fp); 
fclose(fp);


if( (fp=fopen(argv[3], "r")) == NULL ){
    printf("\nError, could not open input file.\n");    
    exit(1);    
}
fread(array2d3, sizeof(array2d3), 1, fp); 
fclose(fp);

if( (fp=fopen(argv[4], "r")) == NULL ){
    printf("\nError, could not open input file.\n");    
    exit(1);    
}
fread(array2d4, sizeof(array2d4), 1, fp); 
fclose(fp);

for (i=0; i<Y_DIM; i++) {
    for (j=0; j<X_DIM; j++) {
        array2d5[i][j] = array2d1[i][j] + scale*array2d2[i][j] - scale*array2d3[i][j] - scale*array2d4[i][j];
    }
}    

printf("writing to output.spn\n");
if( (fp=fopen("output.spn", "w")) == NULL ){
    printf("\nError, could not open output file.\n");    
    exit(1);    
}    
fwrite(array2d5, sizeof(array2d5), 1, fp); 
fclose(fp);

return 0;
}

