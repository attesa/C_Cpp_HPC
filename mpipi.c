//This progtheram is going to calculate pi utlizing MPI
#include "stdio.h"
#include "mpi.h"
#include "unistd.h"
#include "math.h"

int main(int argc, char **argv)
{
  FILE *fp;
  fp = fopen("mpi_1cores.txt","w");
  int rank,size,i,ierr,j,nj,k;
  double pisum[1]={0.0};
  double pires[1]={0.0};
  double tstart,tend;
  int N=1000000;          //rank of summation
  nj=300;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
for(k=0;k<5;k++)
{
for(j=0;j<nj;j++)
{
  pisum[0]=0.0;
  pires[0]=0.0;
  tstart=MPI_Wtime();
  for (i=rank*N/size;i<(rank+1)*N/size;i++) {
    pisum[0]=pisum[0]+pow(-1.0,i)/(2.0*i+1.0);
  }
  MPI_Reduce(pisum,pires,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}
  tend=MPI_Wtime();
  if (rank==0)
  {
    fprintf(fp,"the result is %f \n",4*pires[0]);
    fprintf(fp,"%f , %d ,  %d \n",(tend-tstart),nj,N);
  }
  N=N/10;
  nj=nj*10;
  
}
  ierr=MPI_Finalize();
  return 0;
}  		


  

