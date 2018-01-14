/* This program is to solve the Laplace function in serial fashion */
#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "time.h"

int main (int argc, char *argv[]) {
    time_t start, stop;	
    clock_t begin,end;
    double time_spent;
    double phi_new[40000];
    double phi_old[40000];            //old and new phi 2x2 array
    int i,j,k,l;                //iteration numbers
    int N=80;                      //number of the grids
    int n2;
    int iter;
    double h;                      // step size
    double x,y;
    double diff = 1.0;             //difference 
    double accuracy = 0.00001;
    FILE *fp;
    fp = fopen("result.txt","w");
    h = 1.0/N;  
//    phi_new=(double*)malloc(N*N);
//    phi_old=(double*)malloc(N*N);  //define N*N array, [N*i+j] to use

    for (i=0;i<6;i++) {         //loop for different N 
//        phi_new = calloc(N*N,sizeof(double));
  //      phi_old = calloc(N*N,sizeof(double));  //define N*N array, [N*i+j] to use
	begin = clock();
        h = 1.0/N;
        time(&start);           //time count                           
	for (j=0;j<N;j++) {
	    for (k=0;k<N;k++) {
		phi_new[N*j+k] = 0.0;
		phi_old[N*j+k] = 0.0;             // initial value
	    }          	      
	}
	
	for (j=0;j<N;j++) {
	    phi_old[N*j] = 0.0;
	    phi_old[N*j+N-1] = 0.0;
	}
	for (k=0;k<N;k++) {
	    x = h*k;
	    phi_old[k] =sin(M_PI*x);
	    phi_old[N*(N-1)+k] = sin(M_PI*x)*exp(-M_PI);   //BC 			
	}
        iter = 0;
 	diff = 1.0;
        while (diff > accuracy) {
	    diff = 0.0;
	    for (j=1;j<(N-1);j++) { 
		for (k=1;k<(N-1);k++) {
		    phi_new[N*j+k]=(phi_old[N*(j-1)+k]+phi_old[N*(j+1)+k]+phi_old[N*j+k+1]+phi_old[N*j+k-1])/4.0;
		    diff = diff + pow((phi_new[N*j+k]-phi_old[N*j+k]),2.0); 
		}
	    }                       //iteration
	    
            diff = pow(diff,0.5);	
            
            for (j=1;j<(N-1);j++) {
                for (k=1;k<(N-1);k++) {
                phi_old[N*j+k] = phi_new[N*j+k];
		}
            }

   	    iter+=1;
	}
	end = clock();
	time(&stop);
	time_spent = (double)(end-begin)/CLOCKS_PER_SEC;
 	fprintf(fp,"%d \t %d \t  %f \n",N,iter,time_spent);
	N = N + 20;
//	free(phi_new);
//	free(phi_new);
//	n2=N*N;
//	phi_new=(double*)realloc(phi_new,n2);
//	phi_old=(double*)realloc(phi_old,n2);
   }          
        
}
  
