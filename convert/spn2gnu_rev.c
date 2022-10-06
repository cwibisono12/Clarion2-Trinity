//Matrix Converter into Gnuscope Format
//Useful for analysis of gamma-gamma matrix for Projections


#include<stdio.h>
#include<stdint.h>
#include<stdlib.h>

int main(int argc, char *argv[]) {
printf("starting spn2gnu\n");
//fflush(stdout);
printf("%i\n", argc);
if (argc < 2) {
printf("usage spn2gnu [infile] [outfile]\n");
return -1;
}

FILE *infile = fopen(argv[1], "r"); //this is the spn file
int *matrix = malloc(5000*5000*sizeof(int));
//float *fmatrix;
float *fmatrix = malloc(5000*5000*sizeof(float));
printf("reading matrix\n");

fread(matrix, sizeof(int), 5000*5000, infile);

printf("converting to floats\n");
for (int i=0; i<5000; ++i) {
for (int j=0; j<5000; ++j) {
//printf("%i    %i   \n", i, j);
*(fmatrix + i*5000 + j) = (float)*(matrix + i*5000 + j);
}
}

printf("writing header\n");
//write the header
FILE *outfile = fopen(argv[2], "w");
int four = 12;
int tempbuffer[4];
tempbuffer[0]=5000;
tempbuffer[1]=0;
tempbuffer[2]=0;
tempbuffer[3]=0; 
//int tempbuffer=4096;
int dbuffer=sizeof(int)*5000*5000;
fwrite(&four, sizeof(int), 1, outfile);
fwrite(tempbuffer, 4, 3, outfile);
fwrite(&four, sizeof(int), 1, outfile);
fwrite(&dbuffer, sizeof(int), 1, outfile);
printf("writing data\n");
fwrite(fmatrix, sizeof(float), 5000*5000, outfile);

printf("writing footer\n");
fwrite(&dbuffer, sizeof(int), 1, outfile);

free(matrix);
free(fmatrix);
fclose(outfile);
fclose(infile);

return 0;

}
