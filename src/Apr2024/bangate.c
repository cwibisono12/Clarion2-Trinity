// C function for applying 2D cuts
//J. M Allmond
//data x   data y   ban x array   ban y array   num ban points
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define TRUE 1
#define FALSE 0
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
 
 
 




