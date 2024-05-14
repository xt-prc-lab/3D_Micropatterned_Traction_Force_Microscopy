/*******************************************************************************************************************************
*                                                                                                                              *
*                               Function used to apply a multiplicative 3D window to a 3D stack.                               *
*                                                                                                                              *
*   Inputs:                                                                                                                    *
*       im [3D matrix]: 3D matrix to be windowed.                                                                              *
*       method [int, optional]: type of window applied. 1 means Hanning window, otherwise means no windowing.                  *
*                                                                                                                              *
*   Outputs:                                                                                                                   *
*       im1 [3D matrix]: windowed 3D matrix.                                                                                   *
*                                                                                                                              *
*   Last Revison Date: 04/03/2024                                                                                              *
*   Author: Manuel Gomez Gonzalez                                                                                              *
*                                                                                                                              *
*   References:                                                                                                                *
*       N/A.                                                                                                                   *
*                                                                                                                              *
*******************************************************************************************************************************/

/*-----------------------------------------------------------------------------------------------------------------------------*
*                                                           Imports.                                                           *
*-----------------------------------------------------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <string.h>         /* Needed for memcpy() */
#include "mex.h"
#include "matrix.h"

/*----------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------*
*                               Gateway routine used to interface the C code and the Matlab code.                              *
*-----------------------------------------------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* Declare variables. */
    int method;
    mwSize num_dims;
    const mwSize *dims;
    double *im;
    double *im1;

    /* Check for proper number of input and output arguments. */
    if (nrhs == 1){
        method = 1;
    }
    else if (nrhs == 2){
        method = mxGetScalar(prhs[1]);
    }
    else{
        mexErrMsgTxt("Exactly one or two input arguments required.");
    }

    if (nlhs > 1){
        mexErrMsgTxt("Too many output arguments.");
    }

    /* Check data type of input argument. */
    if (!(mxIsDouble(prhs[0]))){
        mexErrMsgTxt("Input array must be of type double.");
    }

    /* Get the number of elements in the input argument. */
    num_dims = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);

    /* Get the data. */
    im = mxGetPr(prhs[0]);

    /* Get the number of dimensions in the input argument.
    Allocate the space for the return argument */
    plhs[0] = mxCreateNumericArray(num_dims, dims, mxDOUBLE_CLASS, mxREAL);
    im1 = mxGetPr(plhs[0]);

    /* Call the C subroutine. */
    xwindow3D(im, im1, dims[0], dims[1], dims[2], method);

}

/*----------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------*
*                                                                                                                              *
*                              C Function used to apply a multiplicative 3D window to a 3D stack.                              *
*                                                                                                                              *
*   Inputs:                                                                                                                    *
*       im [double*]: pointer to the 3D matrix to be windowed.                                                                 *
*       im1 [double*]: pointer to the windowed 3D matrix.                                                                      *
*       si [unsigned int]: number of rows of the 3D matrices.                                                                  *
*       sj [unsigned int]: number of columns of the 3D matrices.                                                               *
*       sk [unsigned int]: number of layers of the 3D matrices.                                                                *
*       method [unsigned int]: type of window applied. 1 means Hanning window, otherwise means no windowing.                   *
*                                                                                                                              *
*   Outputs:                                                                                                                   *
*       N/A.                                                                                                                   *
*                                                                                                                              *
*   Last Revison Date: 27/03/2024                                                                                              *
*   Author: Manuel Gomez Gonzalez                                                                                              *
*                                                                                                                              *
*   References:                                                                                                                *
*       N/A.                                                                                                                   *
*                                                                                                                              *
*-----------------------------------------------------------------------------------------------------------------------------*/

__inline void xwindow3D(double* im, double* im1, unsigned int si, unsigned int sj, unsigned int sk, unsigned int method) {

    register const double tau_i = 2.*M_PI/((double) si);
    register const double tau_j = 2.*M_PI/((double) sj);
    register const double tau_k = 2.*M_PI/((double) sk);

    register double cos_jk;
    register double cos_k;

    register unsigned int i, j, k;

    double *cos_j;
    double *cos_i;

    if (method == 1){
        
        cos_j = (double *) malloc(sj*sizeof(double));
        cos_i = (double *) malloc(si*sizeof(double));

        for(j=sj; j!=0; j--){
            cos_j[j-1] = 1 - cos(tau_j*((double) j));
        }

        for(i=si; i!=0; i--){
            cos_i[i-1] = 1 - cos(tau_i*((double) i));
        }

        for(k=sk; k!=0; k--){

            cos_k = (1 - cos(tau_k*((double) k)))/8;

            for(j=sj; j!=0; j--){

				cos_jk = cos_j[j-1]*cos_k;

                for(i=si; i!=0; i--){

                    im1[si*sj*(k-1) + si*(j-1) + i-1] = im[si*sj*(k-1) + si*(j-1) + i-1] * cos_i[i-1] * cos_jk;

                }

            }

        }

        free(cos_i);
        free(cos_j);

    }
    else {
        memcpy(im1, im, si*sj*sk*sizeof(double));
    }

    return;

}

/*----------------------------------------------------------------------------------------------------------------------------*/

/******************************************************************************************************************************/