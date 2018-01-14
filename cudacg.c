#include "stdio.h"
#include "unistd.h"
#include "math.h"
#include "stdlib.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
#define N 500
/* cuda version CG solver*/


int main(int argc, char **argv) {
  /* generate matrix and vector first */
  float* A;
  A = (float*)malloc(N*N*sizeof(float));
  float* r;
  float* b;
  float* p;
  float* ap;
  float* x;
  float* ax;
  r = (float*)malloc(N*sizeof(float));
  b = (float*)malloc(N*sizeof(float));
  p = (float*)malloc(N*sizeof(float));
  ap = (float*)malloc(N*sizeof(float));
  x = (float*)malloc(N*sizeof(float));
  ax = (float*)malloc(N*sizeof(float));

  cudaError_t cudaStat;
  cublasStatus_t stat;
  cublasHandle_t handle;


  int i,j;
  for (i=0; i<N*N; i++) {
      A[i]=0.0;
    
  }
  FILE *fp;
  fp = fopen("cgcuda_result.txt","w");





  for (i=1; i<N-1; i++) {
    A[i*N+i] = 2.0;
    A[(i+1)*N+i] =-1.0;
    A[i*N+i+1] =-1.0;

  }
  A[1] = -1.0;
  A[N] = -1.0 ;
  A[0] = 2.0;
  A[(N-1)*N+N-1] = 2.0;
  
  double alpha,beta,rtt,rtt1;
  

  for (i=0;i<N;i++) {
    x[i]= 0.0;//rand()%100/100.0;
    b[i] = i; //(rand() % 1000)/100.0;
    ax[i] = 0.0;
  }

  
  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      ax[i] += A[i*N+j]*x[j];
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
  
  //on the device
  float* d_A;
  float* d_ap;
  float* d_p;
  //
  



  while(rtt > 0.00001) {
  printf("%f rrt %f %f \n", rtt,alpha,beta);

  for (i=0;i<N;i++) {
    ap[i]=0.0;
  }
  
  cudaStat = cudaMalloc((void**)&d_A,N*N*sizeof(*A));
  
  cudaStat = cudaMalloc((void**)&d_ap,N*sizeof(*ap));

  cudaStat = cudaMalloc((void**)&d_p,N*sizeof(*p));

  stat = cublasCreate(&handle);
  stat = cublasSetMatrix(N,N,sizeof(*A),A,N,d_A,N);
  stat = cublasSetVector(N,sizeof(*ap),ap,1,d_ap,1);
  stat = cublasSetVector(N,sizeof(*p),p,1,d_p,1);

  float al = 1.0f;
  float bet = 0.0f;

  stat = cublasSgemv(handle,CUBLAS_OP_N,N,N,&al,d_A,N,d_p,1,&bet,d_ap,1);

  stat = cublasGetVector(N,sizeof(*ap),d_ap,1,ap,1);


  //for (i=0;i<N;i++) {
  //  for (j=0;j<N;j++) {
  //    ap[i] += A[i*N+j]*p[j];
  //  }
  //}
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
  cudaFree(d_A);
  cudaFree(d_p);
  cudaFree(d_ap);
  cublasDestroy ( handle );
  free(A);
  free(r);
  free(b);
  free(p);
  free(ap);
  free(x);
  free(ax);

  
} 
