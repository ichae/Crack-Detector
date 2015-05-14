// (c) Copyright Bing Li 2005 - 2009.
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // L = AM_laplacian(U) where U is 2-D or 3-D matrix    
    
    int i, j, k, idx, ndim, nele;
    const int  *dims;
    double *U, *L;
    
    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    nele = mxGetNumberOfElements(prhs[0]);
    
     /* Check for proper number of input and output arguments. */
    if (nrhs != 1)
        mexErrMsgTxt("One input argument required!");
    if (nlhs != 1)
        mexErrMsgTxt("One ouput argument required!");
    if (ndim>3)
        mexErrMsgTxt("Input matrix must be 2-D or 3-D matrix!");
    
    plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
    
    U = mxGetPr(prhs[0]);    
    L = mxGetPr(plhs[0]);
    
    if (ndim==2)
    {        
        ///////////////////////////////////
        // center square
        ///////////////////////////////////
        for (j=1; j<dims[1]-1; j++)
            for (i=1; i<dims[0]-1; i++)
            {
                idx = j*dims[0]+i;
                L[idx] = (U[idx-1]+U[idx+1]+U[idx+dims[0]]+U[idx-dims[0]])/4 - U[idx];
            }
        ///////////////////////////////////
        // four edges
        ///////////////////////////////////            
        for (i=1; i<dims[0]-1; i++)
        {
            L[i]    = (U[i-1]+U[i+1]+2*U[i+dims[0]])/4 - U[i];
            idx = (dims[1]-1)*dims[0]+i;            
            L[idx]  = (U[idx-1]+U[idx+1]+2*U[idx-dims[0]])/4 - U[idx];
        }
        
        for (j=1; j<dims[1]-1; j++)
        {
            idx = j*dims[0];                        
            L[idx] = (2*U[idx+1]+U[idx+dims[0]]+U[idx-dims[0]])/4 - U[idx];
            idx = (j+1)*dims[0]-1;                        
            L[idx] = (2*U[idx-1]+U[idx+dims[0]]+U[idx-dims[0]])/4 - U[idx];
        }
        ///////////////////////////////////
        // four vertices
        ///////////////////////////////////
        L[0] = (U[1]+U[dims[0]])/2 - U[0];
        L[dims[0]-1] = (U[dims[0]-2]+U[2*dims[0]-1])/2 - U[dims[0]-1];
        L[nele-dims[0]] = (U[nele-dims[0]+1]+U[nele-2*dims[0]])/2 - U[nele-dims[0]];
        L[nele-1] = (U[nele-2]+U[nele-1-dims[0]])/2 - U[nele-1];
    }

    else if (ndim==3)
    {
        int temp = dims[0]*dims[1];
        
        ///////////////////////////////////
        // center cube
        ///////////////////////////////////
        for (k=1; k<dims[2]-1; k++)
            for (j=1; j<dims[1]-1; j++)
                for (i=1; i<dims[0]-1; i++)
                {
                    idx = (k*dims[1]+j)*dims[0]+i;
                    L[idx] = (U[idx-1]+U[idx+1]+U[idx+dims[0]]+U[idx-dims[0]]+U[idx-temp]+U[idx+temp])/6 - U[idx];
                }
        
       
        ///////////////////////////////////
        // six faces
        ///////////////////////////////////
        // k=0 k=end
        for (j=1; j<dims[1]-1; j++)
            for (i=1; i<dims[0]-1; i++)
            {
                idx = j*dims[0]+i;                
                L[idx] = (U[idx-1]+U[idx+1]+U[idx+dims[0]]+U[idx-dims[0]]+2*U[idx+temp])/6 - U[idx];
                idx = ((dims[2]-1)*dims[1]+j)*dims[0]+i;         
                L[idx] = (U[idx-1]+U[idx+1]+U[idx+dims[0]]+U[idx-dims[0]]+2*U[idx-temp])/6 - U[idx];
            }
        // j=0 j=end
        for (j=1; j<dims[2]-1; j++)
            for (i=1; i<dims[0]-1; i++)
            {        
                idx = (j*dims[1])*dims[0]+i;
                L[idx] = (U[idx-1]+U[idx+1]+2*U[idx+dims[0]]+U[idx-temp]+U[idx+temp])/6 - U[idx];
                idx = (j*dims[1]+dims[1]-1)*dims[0]+i;
                L[idx] = (U[idx-1]+U[idx+1]+2*U[idx-dims[0]]+U[idx-temp]+U[idx+temp])/6 - U[idx];
            }
        // i=0 i=end
        for (j=1; j<dims[2]-1; j++)
            for (i=1; i<dims[1]-1; i++)
            {
                idx = (j*dims[1]+i)*dims[0];
                L[idx] = (2*U[idx+1]+U[idx+dims[0]]+U[idx-dims[0]]+U[idx-temp]+U[idx+temp])/6 - U[idx];
                idx = (j*dims[1]+i)*dims[0]+dims[0]-1;
                L[idx] = (2*U[idx-1]+U[idx+dims[0]]+U[idx-dims[0]]+U[idx-temp]+U[idx+temp])/6 - U[idx];
            }
        
        
        ///////////////////////////////////
        // twelve edges
        ///////////////////////////////////
        // i,j = 0,end
        for (i=1; i<dims[2]-1; i++)
        {idx = (k*dims[1]+j)*dims[0]+i;
            idx = (i*dims[1])*dims[0]; 
            L[idx] = (U[idx+1]+U[idx+dims[0]])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
            idx = (i*dims[1])*dims[0]+dims[0]-1; 
            L[idx] = (U[idx-1]+U[idx+dims[0]])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
            idx = (i*dims[1]+dims[1]-1)*dims[0]+dims[0]-1; 
            L[idx] = (U[idx-1]+U[idx-dims[0]])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
            idx = (i*dims[1]+dims[1]-1)*dims[0]; 
            L[idx] = (U[idx+1]+U[idx-dims[0]])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
        }
        // i,k = 0,end
        for (i=1; i<dims[1]-1; i++)
        {
            idx = i*dims[0];
            L[idx] = (U[idx+1]+U[idx+temp])/3+(U[idx-dims[0]]+U[idx+dims[0]])/6 - U[idx];
            idx = i*dims[0]+dims[0]-1;
            L[idx] = (U[idx-1]+U[idx+temp])/3+(U[idx-dims[0]]+U[idx+dims[0]])/6 - U[idx];
            idx = ((dims[2]-1)*dims[1]+i)*dims[0]+dims[0]-1;
            L[idx] = (U[idx-1]+U[idx-temp])/3+(U[idx-dims[0]]+U[idx+dims[0]])/6 - U[idx];
            idx = ((dims[2]-1)*dims[1]+i)*dims[0];
            L[idx] = (U[idx+1]+U[idx-temp])/3+(U[idx-dims[0]]+U[idx+dims[0]])/6 - U[idx];
        }
        // j,k = 0,end
        for (i=1; i<dims[0]-1; i++)
        {
            idx = i;
            L[idx] = (U[idx+dims[0]]+U[idx+temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
            idx = (dims[1]-1)*dims[0]+i;
            L[idx] = (U[idx-dims[0]]+U[idx+temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
            idx = ((dims[2]-1)*dims[1]+dims[1]-1)*dims[0]+i;
            L[idx] = (U[idx-dims[0]]+U[idx-temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
            idx = ((dims[2]-1)*dims[1])*dims[0]+i;
            L[idx] = (U[idx+dims[0]]+U[idx-temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
        }
        ///////////////////////////////////
        // eight vertices
        ///////////////////////////////////        
        L[0] = (U[1]+U[dims[0]]+U[temp])/3 - U[0];
        L[dims[0]-1] = (U[dims[0]-2]+U[2*dims[0]-1]+U[dims[0]-1+temp])/3 - U[dims[0]-1];
        L[dims[0]*(dims[1]-1)] = (U[dims[0]*(dims[1]-1)+1]+U[dims[0]*(dims[1]-2)]+U[dims[0]*(dims[1]-1)+temp])/3 - U[dims[0]*(dims[1]-1)];
        L[temp-1] = (U[temp-2]+U[dims[0]*(dims[1]-1)-1]+U[2*temp-1])/3 - U[temp-1];
        
        L[nele-temp] = (U[nele-temp+1]+U[nele-temp+dims[0]]+U[nele-2*temp])/3 - U[nele-temp];
        L[nele-temp+dims[0]-1] = (U[nele-temp+dims[0]-2]+U[nele-temp+2*dims[0]-1]+U[nele-2*temp+dims[0]-1])/3 - U[nele-temp+dims[0]-1];
        L[nele-dims[0]] = (U[nele-dims[0]+1]+U[nele-2*dims[0]]+U[nele-dims[0]-temp])/3 - U[nele-dims[0]];
        L[nele-1] = (U[nele-2]+U[nele-1-dims[0]]+U[nele-1-temp])/3 - U[nele-1];   
       
    }
}