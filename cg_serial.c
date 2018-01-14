#include "stdio.h"
#include "unistd.h"
#include "math.h"
#include "stdlib.h"
/* CG solver serial version*/

int main(int argc, char **argv) {
  /* generate matrix and vector first */
  int N = 40;
  double A[N][N];
  int i,j;
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      A[i][j]=0.0;
    }
  }
  FILE *fp;
  fp = fopen("cgs.txt","w");





  for (i=1; i<N-1; i++) {
    A[i][i] = 2.0;
    A[i+1][i] =-1.0;
    A[i][i+1] =-1.0;

  }
  A[0][1] = -1.0;
  A[1][0] = -1.0 ;
  A[0][0] = 2.0;
  A[N-1][N-1] = 2.0;
  
  double r[N],b[N],p[N],ap[N],x[N],ax[N];
  double alpha,beta,rtt,rtt1;
  

  for (i=0;i<N;i++) {
    x[i]= 0.0;//rand()%100/100.0;
    b[i] = i; //(rand() % 1000)/100.0;
    ax[i] = 0.0;
  }

  
  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      ax[i] += A[i][j]*x[j];
    }
  }
  
  for (i=0;i<N;i++) {
    r[i] = b[i];
    p[i] = r[i];
  }




  rtt=0.0;
  for (i=0;i<N;i++) {
    rtt += r[i]*r[i];
  }


  while(rtt > 0.00001) {
  printf("%f rrt %f %f \n", rtt,alpha,beta);

  for (i=0;i<N;i++) {
    ap[i]=0.0;
  }

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      ap[i] += A[i][j]*p[j];
    }
  }
  alpha = 0.0;
  for (i=0;i<N;i++) {
    alpha += ap[i]*p[i];
  }
  alpha = rtt/alpha;
  rtt1 = 0.0;
  for (i=0;i<N;i++) {
    x[i] = x[i] + alpha*p[i];
    r[i] = r[i] - alpha*ap[i];
    rtt1 += r[i]*r[i];
  }
  beta = rtt1/rtt;
  rtt = rtt1; 
  for (i=0;i<N;i++) {
    p[i] = r[i] + beta*p[i];
  }
}
  for (i=0;i<N;i++) {
    fprintf(fp, "%f %d \n",x[i],i);
  }

  
} 
