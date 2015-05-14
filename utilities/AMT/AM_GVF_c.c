// (c) Copyright Bing Li 2005 - 2009.
#include "math.h"
#include "mex.h"
#include "matrix.h"

void Lap2(float* U, int Height, int Width, float* L)    // 2D Laplacian
{    
    int i, j, idx, nele = Width*Height;
    for (j=1; j<Width-1; j++)
        for (i=1; i<Height-1; i++)
        {
        idx = j*Height+i;
        L[idx] = (U[idx-1]+U[idx+1]+U[idx+Height]+U[idx-Height])/4 - U[idx];
        }
    
    for (i=1; i<Height-1; i++)
    {
        idx = (Width-1)*Height+i;
        L[i]    = (U[i-1]+U[i+1]+2*U[i+Height])/4 - U[i];
        L[idx]  = (U[idx-1]+U[idx+1]+2*U[idx-Height])/4 - U[idx];
    }
    
    for (j=1; j<Width-1; j++)
    {
        idx = j*Height;
        L[idx] = (2*U[idx+1]+U[idx+Height]+U[idx-Height])/4 - U[idx];
        idx = (j+1)*Height-1;
        L[idx] = (2*U[idx-1]+U[idx+Height]+U[idx-Height])/4 - U[idx];
    }
    
    L[0] = (U[1]+U[Height])/2 - U[0];
    L[Height-1] = (U[Height-2]+U[2*Height-1])/2 - U[Height-1];
    L[nele-Height] = (U[nele-Height+1]+U[nele-2*Height])/2 - U[nele-Height];
    L[nele-1] = (U[nele-2]+U[nele-1-Height])/2 - U[nele-1];
}

void Lap3(float* U, int Height, int Width, int Thick, float* L) // 3D Laplacian
{   
    int temp = Height*Width;
    int nele = temp*Thick;
    int idx, i, j, k;
    ///////////////////////////////////
    // center cube
    ///////////////////////////////////
    for (k=1; k<Thick-1; k++)
        for (j=1; j<Width-1; j++)
            for (i=1; i<Height-1; i++)
    {
        idx = (k*Width+j)*Height+i;
        L[idx] = (U[idx-1]+U[idx+1]+U[idx+Height]+U[idx-Height]+U[idx-temp]+U[idx+temp])/6 - U[idx];
            }
    
    
    ///////////////////////////////////
    // six faces
    ///////////////////////////////////
    // k=0 k=end
    for (j=1; j<Width-1; j++)
        for (i=1; i<Height-1; i++)
    {
        idx = j*Height+i;
        L[idx] = (U[idx-1]+U[idx+1]+U[idx+Height]+U[idx-Height]+2*U[idx+temp])/6 - U[idx];
        idx = ((Thick-1)*Width+j)*Height+i;
        L[idx] = (U[idx-1]+U[idx+1]+U[idx+Height]+U[idx-Height]+2*U[idx-temp])/6 - U[idx];
        }
    // j=0 j=end
    for (j=1; j<Thick-1; j++)
        for (i=1; i<Height-1; i++)
    {
        idx = (j*Width)*Height+i;
        L[idx] = (U[idx-1]+U[idx+1]+2*U[idx+Height]+U[idx-temp]+U[idx+temp])/6 - U[idx];
        idx = (j*Width+Width-1)*Height+i;
        L[idx] = (U[idx-1]+U[idx+1]+2*U[idx-Height]+U[idx-temp]+U[idx+temp])/6 - U[idx];
        }
    // i=0 i=end
    for (j=1; j<Thick-1; j++)
        for (i=1; i<Width-1; i++)
    {
        idx = (j*Width+i)*Height;
        L[idx] = (2*U[idx+1]+U[idx+Height]+U[idx-Height]+U[idx-temp]+U[idx+temp])/6 - U[idx];
        idx = (j*Width+i)*Height+Height-1;
        L[idx] = (2*U[idx-1]+U[idx+Height]+U[idx-Height]+U[idx-temp]+U[idx+temp])/6 - U[idx];
        }
    
    
    ///////////////////////////////////
    // twelve edges
    ///////////////////////////////////
    // i,j = 0,end
    for (i=1; i<Thick-1; i++)
    {
        idx = (i*Width)*Height;
        L[idx] = (U[idx+1]+U[idx+Height])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
        idx = (i*Width)*Height+Height-1;
        L[idx] = (U[idx-1]+U[idx+Height])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
        idx = (i*Width+Width-1)*Height+Height-1;
        L[idx] = (U[idx-1]+U[idx-Height])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
        idx = (i*Width+Width-1)*Height;
        L[idx] = (U[idx+1]+U[idx-Height])/3+(U[idx-temp]+U[idx+temp])/6 - U[idx];
    }
    // i,k = 0,end
    for (i=1; i<Width-1; i++)
    {
        idx = i*Height;
        L[idx] = (U[idx+1]+U[idx+temp])/3+(U[idx-Height]+U[idx+Height])/6 - U[idx];
        idx = i*Height+Height-1;
        L[idx] = (U[idx-1]+U[idx+temp])/3+(U[idx-Height]+U[idx+Height])/6 - U[idx];
        idx = ((Thick-1)*Width+i)*Height+Height-1;
        L[idx] = (U[idx-1]+U[idx-temp])/3+(U[idx-Height]+U[idx+Height])/6 - U[idx];
        idx = ((Thick-1)*Width+i)*Height;
        L[idx] = (U[idx+1]+U[idx-temp])/3+(U[idx-Height]+U[idx+Height])/6 - U[idx];
    }
    // j,k = 0,end
    for (i=1; i<Height-1; i++)
    {
        idx = i;
        L[idx] = (U[idx+Height]+U[idx+temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
        idx = (Width-1)*Height+i;
        L[idx] = (U[idx-Height]+U[idx+temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
        idx = ((Thick-1)*Width+Width-1)*Height+i;
        L[idx] = (U[idx-Height]+U[idx-temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
        idx = ((Thick-1)*Width)*Height+i;
        L[idx] = (U[idx+Height]+U[idx-temp])/3+(U[idx-1]+U[idx+1])/6 - U[idx];
    }
    ///////////////////////////////////
    // eight vertices
    ///////////////////////////////////
    L[0] = (U[1]+U[Height]+U[temp])/3 - U[0];
    L[Height-1] = (U[Height-2]+U[2*Height-1]+U[Height-1+temp])/3 - U[Height-1];
    L[Height*(Width-1)] = (U[Height*(Width-1)+1]+U[Height*(Width-2)]+U[Height*(Width-1)+temp])/3 - U[Height*(Width-1)];
    L[temp-1] = (U[temp-2]+U[Height*(Width-1)-1]+U[2*temp-1])/3 - U[temp-1];
    
    L[nele-temp] = (U[nele-temp+1]+U[nele-temp+Height]+U[nele-2*temp])/3 - U[nele-temp];
    L[nele-temp+Height-1] = (U[nele-temp+Height-2]+U[nele-temp+2*Height-1]+U[nele-2*temp+Height-1])/3 - U[nele-temp+Height-1];
    L[nele-Height] = (U[nele-Height+1]+U[nele-2*Height]+U[nele-Height-temp])/3 - U[nele-Height];
    L[nele-1] = (U[nele-2]+U[nele-1-Height]+U[nele-1-temp])/3 - U[nele-1];
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Fext = AM_GVF(f, mu, ITER, normalize)
    
    int i, d, m, ndim, nele, ITER, normalize;
    int  *dims, Fdims[4];    
    double mu, fmin, fmax, fmm;
    mxArray *rhs[3], *lhs[3], *mxf, *tmp;
    float *f, *fx, *fy, *fz, *u, *G, *lapu, *Fext, mag;
    
    
    ndim = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    nele = mxGetNumberOfElements(prhs[0]);
    
    // Check for proper number of input and output arguments.
    if ((nrhs < 3) || (nrhs > 4))
        mexErrMsgTxt("Invalid input arguments!");
    if (nlhs != 1)
        mexErrMsgTxt("Invalid output arguments!");
    if (ndim>3)
        mexErrMsgTxt("Input matrix must be 2-D or 3-D matrix!");
    
    // inputs
    mxf = mxDuplicateArray(prhs[0]);
    f   = (float *)mxGetData(mxf);
    mu  = *mxGetPr(prhs[1]);
    ITER= *mxGetPr(prhs[2]);
    if (nrhs == 3)
        normalize = 0;
    else
        normalize = *mxGetPr(prhs[3]);
    
    // outputs
    Fdims[0] = dims[0];
    Fdims[1] = dims[1];
    if (ndim==2)        Fdims[2] = 2;
    else if (ndim==3)   { Fdims[2] = dims[2]; Fdims[3] = 3;}
    plhs[0] = mxCreateNumericArray(ndim+1, Fdims, mxSINGLE_CLASS, mxREAL);
    Fext = (float *)mxGetData(plhs[0]);    
  
    // normalize to [0,1]
    fmin = f[0]; fmax = f[0];
    for (i=1;i<nele;i++)
    {
        fmin = (f[i]<fmin)?f[i]:fmin;
        fmax = (f[i]>fmax)?f[i]:fmax;
    }
    fmm = fmax - fmin;
    for (i=0;i<nele;i++)
        f[i] = (f[i]-fmin)/fmm;
    
    /////////////////////////////////////
    // calculate the gradient magnitude square G
    /////////////////////////////////////
   
    G = (float *)mxMalloc(sizeof(float)*nele);
    if (ndim==2)
    {
        mexCallMATLAB(2, lhs, 1, &mxf, "AM_gradient");
        fx = (float *)mxGetData(lhs[0]);
        fy = (float *)mxGetData(lhs[1]);
        for (i=0;i<nele;i++)
        {
            G[i] = fx[i]*fx[i]+fy[i]*fy[i];
        }
    }
    else if (ndim==3)
    {
        mexCallMATLAB(3, lhs, 1, &mxf, "AM_gradient");
        fx = (float *)mxGetData(lhs[0]);
        fy = (float *)mxGetData(lhs[1]);
        fz = (float *)mxGetData(lhs[2]);
        for (i=0;i<nele;i++)
        {
            G[i] = fx[i]*fx[i]+fy[i]*fy[i]+fz[i]*fz[i];
        }
    }
    
    
    /////////////////////////////////////
    // Iteration    
    /////////////////////////////////////
    u = (float *)mxMalloc(sizeof(float)*nele);
    lapu = (float *)mxMalloc(sizeof(float)*nele);
    for (d=0;d<ndim;d++)    // loop over dimensions
    {
        if (d==0) // 1st dimension: u fx   
        {
            for (i=0;i<nele;i++)
                u[i] = fx[i];
            if (ndim==2)
            {
                for (m=0;m<ITER;m++)
                {
                    Lap2(u, dims[0], dims[1], lapu);
                    for (i=0;i<nele;i++)
                        u[i] += mu*lapu[i] - G[i]*(u[i]-fx[i]);
                }
            }
            else
            {
                for (m=0;m<ITER;m++)
                {
                    Lap3(u, dims[0], dims[1], dims[2], lapu);
                    for (i=0;i<nele;i++)
                        u[i] += mu*lapu[i] - G[i]*(u[i]-fx[i]);
                }
            }    
            // format output
            for (i=0;i<nele;i++)            
                Fext[i] = u[i];
        }
        else if (d==1) // 2nd dimension: v fy
        {
            for (i=0;i<nele;i++)
                u[i] = fy[i];
            if (ndim==2)
            {
                for (m=0;m<ITER;m++)
                {
                    Lap2(u, dims[0], dims[1], lapu);
                    for (i=0;i<nele;i++)
                        u[i] += mu*lapu[i] - G[i]*(u[i]-fy[i]);
                }
            }
            else
            {
                for (m=0;m<ITER;m++)
                {
                    Lap3(u, dims[0], dims[1], dims[2], lapu);
                    for (i=0;i<nele;i++)
                        u[i] += mu*lapu[i] - G[i]*(u[i]-fy[i]);
                }
            }    
            // format output
            for (i=0;i<nele;i++)
                Fext[i+nele] = u[i];
        }
        else     // 3rd dimension: w fz
        {
            for (i=0;i<nele;i++)
                u[i] = fz[i];
            {
                for (m=0;m<ITER;m++)
                {
                    Lap3(u, dims[0], dims[1], dims[2], lapu);
                    for (i=0;i<nele;i++)
                        u[i] += mu*lapu[i] - G[i]*(u[i]-fz[i]);
                }
            }    
            for (i=0;i<nele;i++)
                Fext[i+nele*2] = u[i];
        }
    }
        

    /////////////////////////////////////    
    // normalizes the output
    /////////////////////////////////////

    if (normalize)
    {        
        if (ndim==2)
        {
            for (i=0;i<nele;i++)
            {
                mag = sqrt(Fext[i]*Fext[i]+Fext[i+nele]*Fext[i+nele]);
                if (mag)
                {
                    Fext[i]     = Fext[i]/mag;
                    Fext[i+nele]= Fext[i+nele]/mag;
                }
            }
        }
        else
        {
            for (i=0;i<nele;i++)
            {
                mag = sqrt(Fext[i]*Fext[i]+Fext[i+nele]*Fext[i+nele]+Fext[i+nele*2]*Fext[i+nele*2]);
                if (mag)
                {
                    Fext[i]     = Fext[i]/mag;
                    Fext[i+nele]= Fext[i+nele]/mag;
                    Fext[i+nele*2]= Fext[i+nele*2]/mag;
                }
            }
        }
    }
    
    //free the memory
    mxDestroyArray(mxf);
    mxDestroyArray(lhs[0]);
    mxDestroyArray(lhs[1]);    
    mxFree(lapu);
    mxFree(u);
    mxFree(G);            
    if (ndim==3)
    {
        mxDestroyArray(lhs[2]);
    }
    
}