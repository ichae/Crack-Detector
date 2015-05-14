// (c) Copyright Bing Li 2005 - 2009.
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //   [Fx,Fy]   = AM_gradient(F)    when F is 2-D
    // [Fx,Fy,Fz]  = AM_gradient(F)    when F is 3-D
    
    int i, j, k, idx, ndim, nele, temp;
    int  *D, dims[3];
    float *F, *R;
    float *Fx, *Fy, * Fz;
    
    ndim = mxGetNumberOfDimensions(prhs[0]);
    D = mxGetDimensions(prhs[0]);
    dims[0] = D[0];
    dims[1] = D[1];
    if (ndim==2)
        dims[2] = 1;
    else
        dims[2] = D[2];

    temp = dims[0]*dims[1];
    nele = temp*dims[2];

     /* Check for proper number of input and output arguments. */
    if (nrhs > 2)
        mexErrMsgTxt("Invalid input arguments!");
    if (ndim != nlhs)
        mexErrMsgTxt("Invalid output arguments!");
    if (ndim>3)
        mexErrMsgTxt("Input matrix must be 2-D or 3-D matrix!");

    for (i=0; i<ndim; i++)
        plhs[i] = mxCreateNumericArray(ndim, dims, mxSINGLE_CLASS, mxREAL);
    
    F = (float *)mxGetData(prhs[0]);
        
    Fx = (float *)mxGetData(plhs[0]);
    for (k=0; k<dims[2]; k++)
    {
        for (i=0; i<dims[0]; i++)   // boundaries
        {
            Fx[k*temp+i] = 0;
            Fx[(k+1)*temp-1-i] = 0;
        }
        for (j=1; j<dims[1]-1; j++) // interior
            for (i=0; i<dims[0]; i++)
            {
            idx = (k*dims[1]+j)*dims[0]+i;
            Fx[idx] = (F[idx+dims[0]] - F[idx-dims[0]])/2;
            }
    }
    
    Fy = (float *)mxGetData(plhs[1]);
    for (k=0; k<dims[2]; k++)
        for (j=0; j<dims[1]; j++)
        {
            idx = (k*dims[1]+j)*dims[0];
            Fy[idx] = 0;    // boundaries
            for (i=1; i<dims[0]-1; i++) // interior
            {
                idx = (k*dims[1]+j)*dims[0]+i;
                Fy[idx] = (F[idx+1] - F[idx-1])/2;
            }
            idx = (k*dims[1]+j)*dims[0]+dims[0]-1;
            Fy[idx] = 0;    // boundaries
        }

    if (ndim==3)
    {
        Fz = (float *)mxGetData(plhs[2]);
                
        for (j=0; j<dims[1]; j++)   // boundaries
            for (i=0; i<dims[0]; i++)
            {
            idx = j*dims[0]+i;
            Fz[idx] = 0;
            Fz[nele-1-idx] = 0;
            }
        
        for (k=1; k<dims[2]-1; k++) // interior
            for (j=0; j<dims[1]; j++)
                for (i=0; i<dims[0]; i++)
                {
                    idx = (k*dims[1]+j)*dims[0]+i;
                    Fz[idx] = (F[idx+temp] - F[idx-temp])/2;
                }
    }    
    
    if (nrhs == 2)  // anisotropic 3D case
    {
        R = (float *)mxGetData(prhs[1]);
        if (R[0]!=1)
            for(i=0; i<nele; i++)
                Fx[i] = Fx[i]/R[0];
        
        if (R[1]!=1)
            for(i=0; i<nele; i++)
                Fy[i] = Fy[i]/R[1];
                
        if ((R[2]!=1) && (ndim==3))
            for(i=0; i<nele; i++)
                Fz[i] = Fz[i]/R[2];
    }
}
