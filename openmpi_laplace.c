/* This program is to solve the Laplace function with openMP and MPI */
#include "stdio.h"
#include "mpi.h"
#include "unistd.h"
#include "math.h"
#include "stdlib.h"
#include "omp.h"

int main(int argc, char **argv)
{
  
  int rank,size,ierr;
  FILE *fp;
  double local_diff[1]={0.0};
  double diff[1]={0.0};    /* local difference and whole difference*/
  fp = fopen("openmpi_laplace_result_8x50.txt","w");
  MPI_Status status;
  MPI_Request request;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double tstart,tend;
  tstart=MPI_Wtime();  
  /* define 2d array according to number of processor */
  int size_factor = 50;
  int N = size_factor*size;  
  double phi_new[N][N];
  double phi_old[N][N];
  int i,j,isize;
  double h;         /*step size*/
  double tdiff;
  double accuracy = 0.000001;
  int condition[1] = {1};    /* condition for while loop*/
  local_diff[0]=0.0;
  diff[0]=0.0;
  h=1.0/(N-1);
  double x;
  
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      phi_new[i][j]=0.0;
      phi_old[i][j]=0.0;       /* giving initial guess*/	
    }	
  }
  
  for (j=0;j<N;j++) {
    phi_old[0][j] = 0.0;
    phi_old[N-1][j] = 0.0;
  }
  for (i=0;i<N;i++) {
    x = h*i;
    phi_old[i][0] =sin(M_PI*x);
    phi_old[i][N-1] = sin(M_PI*x)*exp(-M_PI);   //BC 			
  } 
  int iter=0;  
  while (condition[0]) {
    /* iteration part part begin*/
    local_diff[0]=0.0;
    diff[0]=0.0;
    tdiff=0.0;
    if (rank==0) 
    {
      #pragma omp parallel for collapse(2) reduction(+:tdiff)
      for (i=1;i<N/2;i++) {
  	for (j=1;j<N-1;j++) {
	  phi_new[i][j] = (phi_old[i-1][j]+phi_old[i+1][j]+phi_old[i][j+1]+phi_old[i][j-1])/4.0;
	  tdiff = tdiff + pow((phi_new[i][j]-phi_old[i][j]),2.0);
	}
      }
      #pragma omp parallel for collapse (2)
      for (i=1;i<N/2;i++) {
        for (j=1;j<N-1;j++) {
          phi_old[i][j]=phi_new[i][j];
        }
      }
      local_diff[0] = tdiff;

    
     
    }   

/*
    else if (rank>0 && rank<size-1) 
    {
      for (i=rank*size_factor;i<(rank+1)*size_factor;i++) {
        for (j=1;j<N-1;j++) {
          phi_new[i][j] = (phi_old[i-1][j]+phi_old[i+1][j]+phi_old[i][j+1]+phi_old[i][j-1])/4.0;
          local_diff[0] = local_diff[0] + pow((phi_new[i][j]-phi_old[i][j]),2.0);
	}
      }
    
      for (i=rank*size_factor;i<(rank+1)*size_factor;i++) {
        for (j=1;j<N-1;j++) {
          phi_old[i][j]=phi_new[i][j];
        }
      }

    }
*/      
    else if(rank==3)
    {
      #pragma omp parallel for collapse(2) reduction (+:tdiff)
      for (i=N/2;i<N-1;i++) {
        for (j=1;j<N-1;j++) {
          phi_new[i][j] = (phi_old[i-1][j]+phi_old[i+1][j]+phi_old[i][j+1]+phi_old[i][j-1])/4.0;
          tdiff = tdiff + pow((phi_new[i][j]-phi_old[i][j]),2.0);
  	}
      }
      #pragma omp parallel for collapse(2) 
      for (i=N/2;i<N-1;i++) {
        for (j=1;j<N-1;j++) {
          phi_old[i][j]=phi_new[i][j];
        }
      }


      local_diff[0] = tdiff;
    }
    /* iteration part over*/
    MPI_Barrier(MPI_COMM_WORLD);
    /* message passing part begin, to avoid blocking, do it sequentially*/
    if (rank==0)
    {
      ierr = MPI_Isend(&phi_old[N/2-1],N,MPI_DOUBLE,3,1,MPI_COMM_WORLD,&request);
      ierr = MPI_Irecv(&phi_old[N/2],N,MPI_DOUBLE,3,1,MPI_COMM_WORLD,&request);
      //ierr = MPI_Isend(phi_old,N-2,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD,&request);
      //ierr = MPI_Irecv(*(*(phi_old+size_factor)+1),N-2,MPI_DOUBLE,rank+1,rank,MPI_COMM_WORLD,&request);
    }
/*    
    else if (rank>0 && rank<size-1)
    {
      ierr = MPI_Irecv(&phi_old[rank*size_factor-1],N,MPI_DOUBLE,rank-1,2*rank-1,MPI_COMM_WORLD,&request);
      ierr = MPI_Isend(&phi_old[rank*size_factor],N,MPI_DOUBLE,rank-1,2*rank-1,MPI_COMM_WORLD,&request);
      //ierr = MPI_Irecv(&phi_old[rank*size_factor-1],N,MPI_DOUBLE,rank-1,rank,MPI_COMM_WORLD,&request);
      ierr = MPI_Isend(&phi_old[(rank+1)*size_factor-1],N,MPI_DOUBLE,rank+1,2*rank+1,MPI_COMM_WORLD,&request);
      ierr = MPI_Irecv(&phi_old[(rank+1)*size_factor],N,MPI_DOUBLE,rank+1,2*rank+1,MPI_COMM_WORLD,&request);
    }   
 for omp */
    else if(rank==3)
    { 
      ierr = MPI_Irecv(&phi_old[N/2-1],N,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&request);
      ierr = MPI_Isend(&phi_old[N/2],N,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&request);

    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* message passing part over */
    
//     if (rank==0) {
//     printf("%f  ",local_diff[0]);
//   }
    /* using allreduce to gather the diff*/
    ierr = MPI_Allreduce(local_diff,diff,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //printf("%f",diff[0]);
    if (pow(diff[0],0.5)<accuracy) {
      condition[0]=0;
    }
    if (rank==0) {
      iter=iter+1;
//      printf("%f  %d \n",diff[0],iter);
    }



  } /* end of while*/

  MPI_Barrier(MPI_COMM_WORLD);
  /* gather data to rank 0 */ 
  int trank; // temporary rank
  for (isize=0;isize<N/2;isize++) {
      if (rank==3) {
	ierr = MPI_Send(&phi_old[isize+N/2],N,MPI_DOUBLE,0,1,MPI_COMM_WORLD);       
      }
      if (rank==0) {
        ierr = MPI_Recv(&phi_old[isize+N/2],N,MPI_DOUBLE,3,1,MPI_COMM_WORLD,&status);
      }
   }
  MPI_Barrier(MPI_COMM_WORLD);

  double phit;    
  tend = MPI_Wtime();
  if (rank == 0) {
    fprintf(fp,"%d  %d %d  %f \n",iter,size,N,tend-tstart);
    for (i=0;i<N;i++) {
      for (j=0;j<N;j++) {
	phit = sin(M_PI*h*i)*exp(-M_PI*h*j);
	fprintf(fp, "%f   %f  %d  %d \n",phi_old[i][j],phit,i,j);

      }
    }
  }

ierr=MPI_Finalize();
return 0;
}
