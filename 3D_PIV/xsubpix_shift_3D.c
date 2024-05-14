/*******************************************************************************************************************************
*                                                                                                                              *
*         Function used to apply a 3D shift to a 3D stack. Subpixel accuracy is achieved through linear interpolation.         *
*                                                                                                                              *
*   Inputs:                                                                                                                    *
*       im [3D matrix]: 3D matrix to be shifted.                                                                               *
*       DX [double]: sub-pixel shift in the x direction, or the columns of the 3D matrix.                                      *
*       DY [double]: sub-pixel shift in the y direction, or the rows of the 3D matrix.                                         *
*       DZ [double]: sub-pixel shift in the z direction, or the layers of the 3D matrix.                                       *
*                                                                                                                              *
*   Outputs:                                                                                                                   *
*       im1 [3D matrix]: shifted 3D matrix.                                                                                    *
*                                                                                                                              *
*   Last Revison Date: 27/03/2024                                                                                              *
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
#include "mex.h"
#include "matrix.h"

/*----------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------*
*                                                            Macros.                                                           *
*-----------------------------------------------------------------------------------------------------------------------------*/

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define sign(a) ( (a > 0) - (a < 0) )

/*----------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------*
*                               Gateway routine used to interface the C code and the Matlab code.                              *
*-----------------------------------------------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Gateway routine used to interface the C code and the Matlab code. */

    /* Declare variables. */
    double Xi ;
    double Xj ;
    double Xk ;
    double padval ;
    mwSize num_dims ;
    const mwSize *dims ;
    double *im ;
    double *im1 ;

    /* Check for proper number of input and output arguments. */
    if (nrhs == 4) {
        padval = 0 ;
    }
    else if (nrhs == 5) {
        padval = mxGetScalar(prhs[4]) ;
    }
    else {
        mexErrMsgTxt("Four or five input argument required.");
    }

    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }

    /* Check data type of input argument. */
    if (!(mxIsDouble(prhs[0]))) {
        mexErrMsgTxt("Input array must be of type double.");
    }

    Xi = mxGetScalar(prhs[2]) ;
    Xj = mxGetScalar(prhs[1]) ;
    Xk = mxGetScalar(prhs[3]) ;

    /* Get the number of elements in the input argument. */
    num_dims = mxGetNumberOfDimensions( prhs[0] );
    dims = mxGetDimensions( prhs[0] ) ;

    /* Get the data. */
    im = mxGetPr( prhs[0] );

    /* Get the number of dimensions in the input argument.
    Allocate the space for the return argument */
    plhs[0] = mxCreateNumericArray( num_dims, dims, mxDOUBLE_CLASS, mxREAL) ;
    im1 = mxGetPr( plhs[0] );

    /* Call the C subroutine. */
    xsubpix_shift_3D( im, im1, dims[0], dims[1], dims[2], Xi, Xj, Xk, padval ) ;

}

/*----------------------------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------------------------------------------*
*                                                                                                                              *
*      C Function used to apply a 3D shift to a 3D stack. Subpixel accuracy is achieved through tri-linear interpolation.      *
*                                                                                                                              *
*   Inputs:                                                                                                                    *
*       im [double*]: pointer to the 3D matrix to be shifted.                                                                  *
*       im1 [double*]: pointer to the shifted 3D matrix.                                                                       *
*       si [unsigned int]: number of rows of the 3D matrices.                                                                  *
*       sj [unsigned int]: number of columns of the 3D matrices.                                                               *
*       sk [unsigned int]: number of layers of the 3D matrices.                                                                *
*       Xi [double]: sub-pixel shift in the y direction, or the rows of the 3D matrix.                                         *
*       Xj [double]: sub-pixel shift in the x direction, or the columns of the 3D matrix.                                      *
*       Xk [double]: sub-pixel shift in the z direction, or the layers of the 3D matrix.                                       *
*       padval [double]: uniform value used in the extrapolation regions.                                                      *
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

__inline void xsubpix_shift_3D( double* im, double* im1, unsigned int si, unsigned int sj, unsigned int sk, 
        double Xi, double Xj, double Xk, double padval ) {

    /* Integer part of the displacement, along the x, y and z directions. */
    register const int I = floorl( -Xi ) ;
    register const int J = floorl( -Xj ) ;
    register const int K = floorl( -Xk ) ;

    const double xi = -Xi - (double) I ;      /* Decimal part of the displacement, i.e. sub-pixel displacement, */
    const double xj = -Xj - (double) J ;      /* along the x, y and z directions. */
    const double xk = -Xk - (double) K ;

    const double fact[8] = { (1.-xi)*(1.-xj)*(1.-xk),
                                 xi *(1.-xj)*(1.-xk),
                             (1.-xi)*    xj *(1.-xk),
                                 xi *    xj *(1.-xk),
                             (1.-xi)*(1.-xj)*    xk ,
                                 xi *(1.-xj)*    xk ,
                             (1.-xi)*    xj *    xk ,
                                 xi *    xj *    xk } ;

    /* Upper and lower bounds of the coordinates between which there is enough information to perform the linear interpolation. */
    register const unsigned int x0 = 1 + min( si, max( 0, -I ) ) ;
    register const unsigned int y0 = 1 + min( sj, max( 0, -J ) ) ;
    register const unsigned int z0 = 1 + min( sk, max( 0, -K ) ) ;
    register const unsigned int xf = max( 0, si + min( 0, -ceill( -Xi ) ) ) ;
    register const unsigned int yf = max( 0, sj + min( 0, -ceill( -Xj ) ) ) ;
    register const unsigned int zf = max( 0, sk + min( 0, -ceill( -Xk ) ) ) ;
    
    register unsigned int i, j, k ;                         /* Counters. */

    register unsigned int n_term = 0 ;                      /* Number of valid terms used in the interpolation. */
    register unsigned int ind ;                             /* Index used to run along fact[ind]. */
    double fact_used[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ;
    unsigned int ind_used[8] = { 0, 0, 0, 0, 0, 0, 0, 0 } ;
    unsigned int term[8] = { 0, 0, 0, 0, 0, 0, 0, 0 } ;

    /* For the coordinate values outisde the x0-xf, y0-yf and z0-zf ranges, there is not enough information to perform a 
     * linear interpolation. For those coordinates, we keep the original values of im. */
    for(k=sk; k>zf; k--){           /* For z > zf. */
        for(j=sj; j>0; j--){
            for(i=si; i>0; i--){
                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = padval ;
            }
        }
    }

    for(k=zf; k>z0-1; k--){

        for(j=sj; j>yf; j--){       /* For y > yf. */
            for(i=si; i>0; i--){
                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = padval ;
            }
        }

        for(j=yf; j>y0-1; j--){

            for(i=si; i>xf; i--){   /* For x > xf. */
                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = padval ;
            }

            /* Tri-linear interpolation. */
            for(i=xf; i>x0-1; i--){

                term[0] = si*sj*(k-1+K) + si*(j-1+J) + (i-1+I) ;
                term[1] = si*sj*(k-1+K) + si*(j-1+J) + (i+I) ;
                term[2] = si*sj*(k-1+K) + si*(j+J) + (i-1+I) ;
                term[3] = si*sj*(k-1+K) + si*(j+J) + (i+I) ;
                term[4] = si*sj*(k+K) + si*(j-1+J) + (i-1+I) ;
                term[5] = si*sj*(k+K) + si*(j-1+J) + (i+I) ;
                term[6] = si*sj*(k+K) + si*(j+J) + (i-1+I) ;
                term[7] = si*sj*(k+K) + si*(j+J) + (i+I) ;

                for(ind=0; ind<7; ind++ ){

                    ind_used[n_term] = ind+1 ;
                    fact_used[n_term] = fact[ ind_used[n_term] ] ;

                    n_term += (int) ( ( term[ind+1] < si*sj*sk )&&( fact_used[n_term]!=0.0 ) ) ;

                }

                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = fact[0] * im[ term[0] ] ;

                for(ind=0; ind<n_term; ind++ ){

                    im1[ si*sj*(k-1) + si*(j-1) + i-1 ] += fact_used[ind] * im[ term[ ind_used[ind] ] ] ;

                }

                n_term = 0 ;

            }

            for(i=x0-1; i>0; i--){  /* For x < x0. */
                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = padval ;
            }

        }

        for(j=y0-1; j>0; j--){      /* For y < y0. */

            for(i=si; i>0; i--){
                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = padval ;
            }

        }

    }

    for(k=z0-1; k>0; k--){          /* For z < z0. */

        for(j=sj; j>0; j--){

            for(i=si; i>0; i--){
                im1[ si*sj*(k-1) + si*(j-1) + i-1 ] = padval ;
            }

        }

    }

    return ;

}

/*----------------------------------------------------------------------------------------------------------------------------*/

/******************************************************************************************************************************/