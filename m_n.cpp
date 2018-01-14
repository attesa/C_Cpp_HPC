#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <mpi.h>
using namespace std;
// point to point implementation of mandelbrot set problem
/*
 * mpicxx m_n.cpp  -std=c++11  -DMPICH_IGNORE_CXX_SEEK
 * mpirun -n 8 ./a.out
 *
 */
int main()
{
	MPI_Status status;
	MPI_Request request;
        double tstart,tend;

        MPI_Init(NULL, NULL);
        int rank,size,ierr;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
	//define dimension grid
	int N = 1000;
	//define flag for convergence 
	int conver[N][N];
	int itermax=10000; //interation max steps
	

	//parallel computing
        tstart=MPI_Wtime();	
	for (int i = 0; i < N; ++i) {
	    for (int j = 0; j < N; ++j) {
		conver[i][j] = 0;
	    }
	}
	
	int i_rank; 
	int row;
	int col;
	for (int i=0; i< N*N; i=i+size) {
	    int i_rank = i+rank;
	    row = i_rank/N;
	    col = i_rank%N;
	    if (rank==0) 
	    {
		cout << rank << "\t" << row << "\t" << col << "\t" <<  endl;
	    } 
            double c_im = (col-N/2.0)*4.0/N; // match the grid into
            double c_re = (row-N/2.0)*4.0/N; // a (-2,2) X (-2,2) complex plane
            double x = 0; //initial point(x,y)
            double y = 0;
            int iter = 0;
            while(x*x+y*y<=4 && iter <itermax){     //iteration steps
                    double x_new = x*x - y*y + c_re;
                    y = 2*x*y + c_im;
                    x = x_new;
                    iter++;
            }
            if (iter == itermax)
            {
                conver[row][col] = 1;
            }
	
	}
	
	//End of computing
	//

        //cout <<  "\t" << rank << "\t" <<  endl;	
        MPI_Barrier(MPI_COMM_WORLD);
	//message passing using reduce
	int final_conver[N][N]; 
	ierr = MPI_Reduce(&conver[0][0],&final_conver[0][0],N*N,MPI_INT,MPI_SUM, 0, MPI_COMM_WORLD);
	tend=MPI_Wtime();
	ierr=MPI_Finalize();

	if (rank == 0)
	{
	    cout<<tend-tstart<<endl;
	    ofstream mdbfile;
            mdbfile.open ("wh_s4.txt");
	    //print all
	    /*
	    for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
                mdbfile << (i-N/2.0)*4.0/N << "\t" << (j-N/2.0)*4.0/N << "\t"<< final_conver[i][j] << "\t" <<  endl;
	        }
	    }
	    */	


            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
		    if (final_conver[i][j] > 0 )
		    {
                    mdbfile << (i-N/2.0)*4.0/N << "\t" << (j-N/2.0)*4.0/N <<   endl;
		    }
                }
            }

	}


	return 0;
}


/*
	 for( int row = 0; row < N; ++row ) {
	 	for (int col = 0; col < N; ++col)
	 	{
	 		double c_re = (col-N/2.0)*4.0/N; // match the grid into 
	 		double c_im = (row-N/2.0)*4.0/N; // a (-2,2) X (-2,2) complex plane
	 		double x = 0; //initial point(x,y)
			double y = 0;
			int iter = 0;
			while(x*x+y*y<=4 && iter <itermax){	//iteration steps
				double x_new = x*x - y*y + c_re;
				y = 2*x*y + c_im;
				x = x_new;
				iter++;
			}
			if (iter < itermax)
			{
			}
			else;
	 	}
      
   }

return (0);
}
*/
