#include "mex.h"
#include <math.h>

#include <iostream>
#include <complex>
#include <vector>
#include <limits>  
#ifdef NAN
/* NAN is supported */
#endif
using namespace std;
void quadsolve(vector<complex<double>>& A, vector<complex<double>>& B, vector<complex<double>>& C, double *outMatrix_x1_real, double *outMatrix_x1_imag, double *outMatrix_x2_real, double *outMatrix_x2_imag, double *outMatrix_D_real, double *outMatrix_D_imag)
{
    
    for (size_t i = 0; i < A.size(); i++)
        {
            double eps =0.0001;
            const complex<double> a = A[i], b = B[i], c = C[i];
           
                if (abs(a)+abs(b) < eps && abs(c)>eps) {
                    outMatrix_x1_real[i] = NAN;
                outMatrix_x1_imag[i] = NAN;
                outMatrix_x2_real[i] = NAN;
                outMatrix_x2_imag[i] = NAN;
                    
                    mexPrintf("no solution in position ");
                    char buff[10];
                    sprintf(buff,"%d",i+1);
                    mexPrintf(buff);
                    mexPrintf("\n");
                    continue;
                    }
                 if (abs(a)+abs(b) < eps && abs(c) < eps) {
                    outMatrix_#include "mex.h"
#include <math.h>

#include <iostream>
#include <complex>
#include <vector>
#include <limits>  
#ifdef NAN
/* NAN is supported */
#endif
using namespace std;
void quadsolve(vector<complex<double>>& A, vector<complex<double>>& B, vector<complex<double>>& C, double *outMatrix_x1_real, double *outMatrix_x1_imag, double *outMatrix_x2_real, double *outMatrix_x2_imag, double *outMatrix_D_real, double *outMatrix_D_imag)
{
    
    for (size_t i = 0; i < A.size(); i++)
        {
            double eps =0.0001;
            const complex<double> a = A[i], b = B[i], c = C[i];
           
                if (abs(a)+abs(b) < eps && abs(c)>eps) {
                    outMatrix_x1_real[i] = NAN;
                outMatrix_x1_imag[i] = NAN;
                outMatrix_x2_real[i] = NAN;
                outMatrix_x2_imag[i] = NAN;
                    
                    mexPrintf("no solution in position ");
                    char buff[10];
                    sprintf(buff,"%d",i+1);
                    mexPrintf(buff);
                    mexPrintf("\n");
                    continue;
                    }
                 if (abs(a)+abs(b) < eps && abs(c) < eps) {
                    outMatrix_x1_real[i] = NAN;
                outMatrix_x1_imag[i] = NAN;
                outMatrix_x2_real[i] = NAN;
                outMatrix_x2_imag[i] = NAN;
                    
                    mexPrintf("all solution in position ");
                    char buff[10];
                    sprintf(buff,"%d",i+1);
                    mexPrintf(buff);
                    mexPrintf("\n");
                    continue;
                    }
            
            const complex<double> d = b * b - double(4) * a * c;
            outMatrix_D_real[i] = d.real();
            outMatrix_D_imag[i] = d.imag();
            complex<double> x1;
            complex<double> x2;
            double epsilon = 1e-5;// for comparing doubles
            if (abs(a) > epsilon)
            {
                x1 = (-b - sqrt(d)) / (double(2) * a);
                x2 = (-b + sqrt(d)) / (double(2) * a);
            }
            else if (abs(b) > epsilon)
            {
                x1 = -c / b;
                x2 = -c / b;
            }
            else
            {
                complex<double> x1 = numeric_limits<double>::quiet_NaN();
                complex<double> x2 = numeric_limits<double>::quiet_NaN();
            }
                if (x1 != x2) {
            outMatrix_x1_real[i] = x1.real();
            outMatrix_x1_imag[i] = x1.imag();
            outMatrix_x2_real[i] = x2.real();
            outMatrix_x2_imag[i] = x2.imag();
                }
            if (x1 == x2) {
            outMatrix_x1_real[i] = x1.real();
            outMatrix_x1_imag[i] = x1.imag();
            outMatrix_x2_real[i] = NAN;
            outMatrix_x2_imag[i] = NAN;
            mexPrintf("one solution in position ");
                    char buff[10];
                    sprintf(buff,"%d",i+1);
                    mexPrintf(buff);
                    mexPrintf("\n");
                }
                
            
        }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[])
{
    
    size_t nrows, ncols;                   /* size of matrices */
    /* input matrices */
    vector<complex<double>> A, B, C;
    /* output matrices */
    double *outMatrix_x1_real, *outMatrix_x1_imag, *outMatrix_x2_real, *outMatrix_x2_imag, *outMatrix_D_real, *outMatrix_D_imag;
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:quadsolve:nrhs","Three inputs required.");
    }
    if(nlhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:quadsolve:nlhs","3 outputs required.");
    }
    if( (mxGetM(prhs[0])!= mxGetM(prhs[1])) || (mxGetM(prhs[1])!= mxGetM(prhs[2]))
    ||(mxGetN(prhs[0])!= mxGetN(prhs[1])) || (mxGetN(prhs[1])!= mxGetN(prhs[2])) ) {
        mexErrMsgIdAndTxt("MyToolbox:quadsolve:nrhs","Input matrices must be the same size.");
    }
 
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
   
    
    double *real_data_ptr_A = (double *)mxGetPr(prhs[0]);
    double *imag_data_ptr_A = (double *)mxGetPi(prhs[0]);
    double *real_data_ptr_B = (double *)mxGetPr(prhs[1]);
    double *imag_data_ptr_B = (double *)mxGetPi(prhs[1]);
    double *real_data_ptr_C = (double *)mxGetPr(prhs[2]);
    double *imag_data_ptr_C = (double *)mxGetPi(prhs[2]);
    
    
    for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j < ncols; j++)
        {
            A.push_back(complex<double>(*real_data_ptr_A, *imag_data_ptr_A));
            B.push_back(complex<double>(*real_data_ptr_B, *imag_data_ptr_B));
            C.push_back(complex<double>(*real_data_ptr_C, *imag_data_ptr_C));
            real_data_ptr_A++;
            imag_data_ptr_A++;
            real_data_ptr_B++;
            imag_data_ptr_B++;
            real_data_ptr_C++;
            imag_data_ptr_C++;
        }
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxCOMPLEX);
    outMatrix_x1_real = mxGetPr(plhs[0]);
    outMatrix_x1_imag = mxGetPi(plhs[0]);
    outMatrix_x2_real = mxGetPr(plhs[1]);
    outMatrix_x2_imag = mxGetPi(plhs[1]);
    outMatrix_D_real = mxGetPr(plhs[2]);
    outMatrix_D_imag = mxGetPi(plhs[2]);
    
    /* call the computational routine */
    quadsolve(A, B, C, outMatrix_x1_real, outMatrix_x1_imag, outMatrix_x2_real, outMatrix_x2_imag, outMatrix_D_real, outMatrix_D_imag);
}