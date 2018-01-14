/*Using conjugate gradient method to solve linear equations with mpi implementation */

#include "stdio.h"
#include "mpi.h"
#include "unistd.h"
#include "math.h"
#include "stdlib.h"
int main(int argc, char **argv) {
  /* generate matrix and vector first */  
  int N = 500;
  double A[N][N];
  FILE *fp;
  fp = fopen("cg.txt","w");
  double b[N];
  /* solve equation Ax=b */
  double r[N];
  double p[N];
  double x[N];
  int i,j,k;
  int scale;
  double r_sum[1];
  double alpha_local_den[1];
  double alpha_total_den[1];
  double alpha[1];
  double r1_local_sum[1];
  double r1_sum[1];
  double beta[1];
  double x_out[N];
  double r_local_sum[1];
  double ap[N];
  double ap_r[N];
  double tstart,tend;
  double mts,mte,mtime;
  /* a simple laplacian matrix */

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      A[i][j]=0.0;
    }
  }




  for (i=1; i<N-1; i++) {
    A[i][i] = 2.0; 
    A[i+1][i] =-1.0; 
    A[i][i+1] =-1.0;
     
  }
  A[0][1] = -1.0;
  A[1][0] = -1.0;
  A[0][0] = 2.0; 
  A[N-1][N-1] = 2.0;
 
  /* initial guess for x is 0 */    
  for (i=0; i<N; i++) {
    x[i] = 0.0;
  }
  /* for first step, r=b=p */
  for (i=0; i<N; i++) {
   // b[i] = 1.0;
    b[i] = i; //(rand() % 1000)/100.0;
  }
  for (i=0; i<N; i++) {
    r[i] = b[i];
    p[i] = b[i];
  } 
  int rank,size,ierr;
  MPI_Status status;
  MPI_Request request;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  mtime=0.0;
  scale = N/size;//decide block size depending on number of processors
  tstart=MPI_Wtime();
  r_local_sum[0] = 0.0;
  for (i = rank*scale; i<(rank+1)*scale; i++) {
    r_local_sum[0] += r[i]*r[i];
  }
  MPI_Allreduce(r_local_sum,r_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  printf("%f \n \n",r_sum[0]);
  //k = 0;

  while(r_sum[0]>0.000001) {
     // k += 1;
    /* calculate A*P */
/*   if (rank == 0) {
      ap[0] = A[0][0]*p[0] + A[0][1]*p[1]; 
      for (i=1; i<scale; i++) {
	ap[i] = A[i][i-1]*p[i-1] + A[i][i]*p[i] + A[i][i+1]*p[i+1]; 
      }       
    }
    else if (rank == size-1) {
      for (i=rank*scale; i<N-1; i++) {
        ap[i] = A[i][i-1]*p[i-1] + A[i][i]*p[i] + A[i][i+1]*p[i+1];
      }
      ap[N-1] = A[N-1][N-1]*p[N-1] + A[N-1][N-2]*p[N-2];
    }
    else {
    
*/  mts = MPI_Wtime();
    for (i=rank*scale; i<(rank+1)*scale; i++) {
      ap[i] = 0;
    }
    for (i=rank*scale; i<(rank+1)*scale; i++) {
      for (j=0;j<N;j++) {
        ap[i] += A[i][j]*p[j];
      }
    } 
    mte = MPI_Wtime();
    mtime += mte-mts;

/*      for (i=rank*scale; i<(rank+1)*scale; i++) {
        ap[i] = A[i][i-1]*p[i-1] + A[i][i]*p[i] + A[i][i+1]*p[i+1];
      }
    }*/
//Above can choose a optimized algorithm for Laplacian or a general 
//Matrix Vector multiplication
    

    /* pass a*p to root */
    MPI_Gather(&ap[rank*scale],scale,MPI_DOUBLE,ap_r,scale,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
   //root deal with all the rest of the iteration, which is trivial comparing to matrix vector multiply
    if (rank == 0) {
//      for (i = 0; i<N; i++) {
//        printf("%f %d\n", ap[i],i);
//      }

    
      alpha_total_den[0] = 0.0;
      /* alpha */
      for (i = 0; i<N; i++) {
        alpha_total_den[0] += p[i]*ap_r[i];
      }
      alpha[0] = r_sum[0] / alpha_total_den[0];
//      printf("%f alpha \n",alpha[0]);     
    
    /* Generate new x,r,p */ 
      for (i = 0; i<N; i++) {
        x[i] = x[i] + alpha[0]*p[i];
        r[i] = r[i] - alpha[0]*ap_r[i];
      }
      r1_sum[0] = 0.0;
      for (i =0 ; i<N; i++) {
        r1_sum[0] += r[i]*r[i];
      }
      beta[0] = r1_sum[0]/r_sum[0];
      r_sum[0] = r1_sum[0];
 
      for (i = 0; i<N; i++) {
        p[i] = r[i] + beta[0]*p[i];
      }
    }
    //send r and r_sum back to all processors
    MPI_Bcast(p,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(r_sum,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
     
  }/* end of iteration */
  tend=MPI_Wtime();

  if (rank == 0 ) {
  fprintf(fp, "%f \n",tend-tstart);
  fprintf(fp, "%f matrix time\n",mtime);
    for (i=0;i<N;i++) {
      fprintf(fp, "%f %d \n",x[i],i);
    } 
  }



  ierr=MPI_Finalize();
  return 0;
} 

 
