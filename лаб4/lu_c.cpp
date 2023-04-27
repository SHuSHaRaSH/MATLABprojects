#include "mex.h"
#include <iostream>
#include <complex>
#include <vector>
#include <limits>  
using namespace std;
void lu_c(vector <double>& A, double *outMatrix_L, 
		double *outMatrix_U, int n)
{
    double eps  = 0.0001;
	for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            outMatrix_U[i + j*n] = 0;
            outMatrix_L[i + j*n] = 0;
        }
    }
	for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            outMatrix_U[i*n] = A[i*n];
            outMatrix_L[i] = A[i] / outMatrix_U[0];
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += outMatrix_L[i + k*n] * outMatrix_U[k + j*n];
            }
            outMatrix_U[i + j*n] = A[i + j*n] - sum;
            if (i > j)
            {
                outMatrix_L[j + i*n] = 0;
            }
            else
            {
                if (abs(outMatrix_U[i*n + i]) < eps)
                    mexErrMsgIdAndTxt("MyToolbox:lu_c:nrhs","There is no LU-decomposition");
                else {
                    sum = 0;
                    for (int k = 0; k < i; k++)
                    {
                        sum += outMatrix_L[j + k*n] * outMatrix_U[k + i*n];
                    }
                    outMatrix_L[j + i*n] = (A[j + i*n] - sum) / outMatrix_U[i*n + i];
                }
            }
        }
    }

}

void proisv(vector <vector <double>> A, vector <vector <double>> B, 
			vector <vector <double>> &R, int n)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			for(int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}

void show(vector <vector <double>> A, int n)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			cout <<"\t"<< A[i][j] << "\t";
		}
		cout << endl;
	}
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[])
{
    size_t nrows, ncols, n;                   /* size of matrices */
    /* input m#include "mex.h"
#include <iostream>
#include <complex>
#include <vector>
#include <limits>  
using namespace std;
void lu_c(vector <double>& A, double *outMatrix_L, 
		double *outMatrix_U, int n)
{
    double eps  = 0.0001;
	for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            outMatrix_U[i + j*n] = 0;
            outMatrix_L[i + j*n] = 0;
        }
    }
	for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            outMatrix_U[i*n] = A[i*n];
            outMatrix_L[i] = A[i] / outMatrix_U[0];
            double sum = 0;
            for (int k = 0; k < i; k++)
            {
                sum += outMatrix_L[i + k*n] * outMatrix_U[k + j*n];
            }
            outMatrix_U[i + j*n] = A[i + j*n] - sum;
            if (i > j)
            {
                outMatrix_L[j + i*n] = 0;
            }
            else
            {
                if (abs(outMatrix_U[i*n + i]) < eps)
                    mexErrMsgIdAndTxt("MyToolbox:lu_c:nrhs","There is no LU-decomposition");
                else {
                    sum = 0;
                    for (int k = 0; k < i; k++)
                    {
                        sum += outMatrix_L[j + k*n] * outMatrix_U[k + i*n];
                    }
                    outMatrix_L[j + i*n] = (A[j + i*n] - sum) / outMatrix_U[i*n + i];
                }
            }
        }
    }

}

void proisv(vector <vector <double>> A, vector <vector <double>> B, 
			vector <vector <double>> &R, int n)
{
	for(int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
			for(int k = 0; k < n; k++)
				R[i][j] += A[i][k] * B[k][j];
}

void show(vector <vector <double>> A, int n)
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < n; j++)
		{
			cout <<"\t"<< A[i][j] << "\t";
		}
		cout << endl;
	}
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[])
{
    size_t nrows, ncols, n;                   /* size of matrices */
    /* input matrices */
    //double * A = new double;
    /* output matrices */
    double *outMatrix_L, *outMatrix_U;
    /* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:lu_c:nrhs","One input required.");
    }
    if(nlhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:lu_c:nlhs","Two outputs required.");
    }
    if( (mxGetM(prhs[0])!= mxGetN(prhs[0])) ) {
        mexErrMsgIdAndTxt("MyToolbox:lu_c:nrhs","Input matrice must be square.");
    }
    
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    n = nrows;
    //double * A = new double [nrows, ncols];
    vector <double> A;
    double *real_data_ptr_A = (double *)mxGetPr(prhs[0]);
    for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j < ncols; j++)
        {
           A.push_back(*real_data_ptr_A);
           //A[i*nrows + j] = real_data_ptr_A[i*nrows + j];
           real_data_ptr_A++;
        }
    
    /* create the output matrices */
    for (size_t i = 0; i < nrows; i++) {
        if (A[0] == 0)
            mexErrMsgIdAndTxt("MyToolbox:lu_c:nrhs","There is no LU-decomposition");
    }
    plhs[0] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxREAL);

    
    outMatrix_L = mxGetPr(plhs[0]);
    outMatrix_U = mxGetPr(plhs[1]);


    
    /* call the computational routine */
    lu_c(A, outMatrix_L, outMatrix_U, n);
}