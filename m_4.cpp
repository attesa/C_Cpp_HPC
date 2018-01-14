#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>
#include <mpi.h>
using namespace std;
//4 regions mpi implementation of mandelbrot set problem 
/*
 *
 * mpicxx m_4.cpp  -std=c++11  -DMPICH_IGNORE_CXX_SEEK
 * mpirun -n 4 ./a.out




 */


int main()
{
// Todo place code here
	MPI_Status status;
	MPI_Request request;
        double tstart,tend;

        MPI_Init(NULL, NULL);
        int rank,size,ierr;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
	//define dimension grid
	int N = 1000;
	double grid_row[N];
	double grid_col[N];
	//define flag for convergence 
	int conver[N][N];
	int itermax=10000; //interation max steps
	
        if (size != 4)
        {
                cout << "there are " << size << "process, need 4 process" << endl;
                return 1;
        }

	//parallel computing
        tstart=MPI_Wtime();	
	for (int i = 0; i < N; ++i) {
	    for (int j = 0; j < N; ++j) {
		conver[i][j] = 0;
	    }
	}


	MPI_Barrier(MPI_COMM_WORLD);
	for (int row = 0; row < N/2; ++row) {
	    for (int col = 0; col < N/2; ++col)
	    {
		//give initial value
                double c_im = (col-(rank%2)*N/2.0)*4.0/N; // match the grid into
                double c_re = (row-(rank/2)*N/2.0)*4.0/N; // a (-2,2) X (-2,2) complex plane
                double x = 0; //initial point(x,y)
                double y = 0;
		int iter = 0;
                while(x*x+y*y<=4 && iter <itermax){     //iteration steps
                        double x_new = x*x - y*y + c_re;
                        y = 2*x*y + c_im;
                        x = x_new;
                        iter++;
                }
                grid_row[row] = c_im;
                grid_col[col] = c_re;

		if (iter == itermax)
		{
		    conver[row+(1-rank/2)*N/2][col+(1-rank%2)*N/2] = 1;
		}
		//else
		//{
	        //    conver[row+(1-rank/2)*N][col+(1-rank%2)*N] = 1;
		//}
		//if (rank == 2) 
		//{
                //cout << rank << "\t" << c_re << "\t" << c_im << "\t" <<  endl;
		//}
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
            mdbfile.open ("wh_4.txt");
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
