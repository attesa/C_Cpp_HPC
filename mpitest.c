/* This program is for ping-pong testing */
#include "stdio.h"
#include "mpi.h"
#include "unistd.h"
int main(int argc, char **argv)
{
  
  FILE *fp;
  fp = fopen("mpi_test_result_interti.txt","w");
  int i,j,ierr,rank,size,dest1,dest2,source1,source2,count,tag1,tag2;
  int stat_count,stat_source,stat_tag;
  double data[1000000];
  double time1,time2,time;
  MPI_Status status;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("I am process %d of %d\n",rank,size);
  if (size!=2) {
    return 1;
  }  
  count=1000;
  source1=0;
  source2=1;
  dest1=0;
  dest2=1;
  tag1=6;
  tag2=7;
for (j=1;j<=4;j++) {
  time1=MPI_Wtime();  
  for (i=1;i<=1000;i++){
    /* send and receive for p1 */
    if (rank==source1) {
     ierr = MPI_Send(data,count,MPI_DOUBLE,dest2,tag1,MPI_COMM_WORLD);
     ierr = MPI_Recv(data,count,MPI_DOUBLE,source2,tag2,MPI_COMM_WORLD,&status);
    }

    else if (rank==source2) {
     ierr = MPI_Recv(data,count,MPI_DOUBLE,source1,tag1,MPI_COMM_WORLD,&status);
     ierr = MPI_Send(data,count,MPI_DOUBLE,dest1,tag2,MPI_COMM_WORLD);
    } 
  } 
  time2=MPI_Wtime();
  printf("%f %f",time1,time2);
  time=(time2-time1)/2000.0;   
  ierr = MPI_Get_count(&status,MPI_DOUBLE,&stat_count);
  printf ("Status of receive : dest= %d, source=%d, tag=%d,\
  count=%d\n",rank,stat_source,stat_tag,stat_count);
  ierr = MPI_Get_count(&status,MPI_DOUBLE,&stat_count);
  fprintf(fp,"%d %f\n",count*8,time);
  count=count*10;
}
  ierr=MPI_Finalize();
  return 0;
}
