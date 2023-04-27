v01.00v00.00 ��  �  �  ��QZ����U�;p|C��9�?Q�3>����3�rV��pa�Au�����K�P"�N@��Rd�I̍k7��42<���1q�,5��m<�fњ�X}���LG���Fi:��S�Ř�;�]�UΞ�� � �dP�شo�TlJ�4S<���5b�q"�>���m'e2)X�T�Sԛ3P C�:|���K�o�o���_Q���w�'�47�V��X�X*��?'B��I�ӃvY� M�(S��=�O�� ]�f4�!ſ6(��В�YoOTh=R�e��>lKÙy�U���K;�ZH|tte�5���K!����l�eY0��)vL,-L��p���������C_^��X߼�s�|z�ޙ�[&2�>�Y����5Ǚ�����u�1��Ig��T@��j��y�`~�o�����D9G�k����!H��y1G�+3�F}�#�jYZY�	��'h`��
RN�Ä��EJ�>hk��!���C����3������'�&��V:}Ï�&W���b���S8�%�J��
v�3R�To�,�s�9~|&3V�$�E<Ɓ���Ӫ���ë�I��GB8��d)h�Y털Ng/��	,�w���4}�w5o��x��N�Ž�ͺ2R}�_ �	��V�J��VC�s����>�+��N��B�eG���2�,��_LJ�8ԷV�>��U�&���7Y�F�ey������9���`o�6T94Mz��L��,Ep��s���a��2��#include "mex.h"
#include <iostream>
#include <complex>
#include <vector>
#include <limits>  
using namespace std;
void quadsolve(vector <double> & A, double *outMatrixL, double *outMatrixU)
{
     for (size_t i = 0; i < A.size(); i++)
        {
	    const double a = A[i];
	    
            outMatrixU[i] = a;
            outMatrixL[i] = a;}
        
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs,  const mxArray *prhs[])
{
    size_t nrows, ncols;                   /* size of matrices */
    	/* input matrices */
    vector <double> A; 
    	/* output matrices */
    double *outMatrixL, *outMatrixU;
    	/* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:quadsolve:nrhs","1 inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:quadsolve:nlhs","2 outputs required.");
    }
   /* if( (mxGetM(prhs[0])!= mxGetM(prhs[1])) || (mxGetM(prhs[1])!= mxGetM(prhs[2]))
    ||(mxGetN(prhs[0])!= mxGetN(prhs[1])) || (mxGetN(prhs[1])!= mxGetN(prhs[2])) ) {
        mexErrMsgIdAndTxt("MyToolbox:quadsolve:nrhs","Input matrices must be the same size.");
    }
    */
    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    

    
    
 /*   for (size_t i = 0; i < nrows; i++)
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
        }*/
    
    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize) nrows, (mwSize) ncols, mxREAL);
   
    
    /* call the computational routine */
    quadsolve(prhs[0], plhs[0], plhs[1]);
}