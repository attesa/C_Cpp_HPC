/* This program is an OpenMP example for calculate pi */
#include "stdio.h"
#include "omp.h"
#include "math.h"

int main (int argc, char *argv[]) {
FILE *fp;
fp = fopen("omp_1cores.txt","w");
int i,n,j,nj,nu,k;
double a;
a=0.0;
n=1000000;  //total summation number
nj=30;      //loop number for loop j
nu=1;
double tstart,tend;
for (k=0;k<5;k++){
omp_set_num_threads(nu);
tstart = omp_get_wtime();
for (j=0;j<nj;j++){
a=0.0;
#pragma omp parallel for private(i) shared(n) \
reduction(+:a) 
for (i=0;i<n;i++) {
a +=pow(-1.0,i)/(2.0*i+1.0) ;
}

}
tend=omp_get_wtime();
fprintf(fp,"pi=%f \n",4*a);
fprintf(fp,"Work took %f seconds for n=%d \n", (tend - tstart)/nj,n);
n=n/10;
nj=nj*10;
}
}
