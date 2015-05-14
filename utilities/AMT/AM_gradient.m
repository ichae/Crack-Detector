function varargout = AM_gradient(A, R)
% AM_GRADIENT    2D/3D gradient, a faster implementation than build-in GRADIENT().
%     [Fx,Fy]     = AM_gradient(F)    when F is a 2-D matrix
%     [Fx,Fy,Fz]  = AM_gradient(F)    when F is a 3-D matrix
%     [...]       = AM_gradient(F, R) when F is an anisotropic data
% 
%     Returns the numerical gradient of the matrix F. Fx corresponds to dF/dx,
%     the differences in the x (column) direction. Fy corresponds to dF/dy,
%     the differences in the y (row) direction. Fz corresponds to dF/dz.
%     Fx, Fy and Fz are the same size as F. R is a N-element vector, [Rx Ry
%     Rz] define the spacing between point in F for x, y and z dimension,
%     default value is [1 1 1].  
% 
%     Note that for memory saving reason, the output data class is
%     single/float (32-bit), and the results at the boundaries are
%     different from the ones generated from GRADIENT().   
% 
%     Example
%     [fx,fy] = AM_gradient(magic(5));
% 
%     fx =
% 
%              0   -8.0000   -8.0000    7.0000         0
%              0   -8.0000    4.5000    4.5000         0
%              0    4.5000    7.0000    4.5000         0
%              0    4.5000    4.5000   -8.0000         0
%              0    7.0000   -8.0000   -8.0000         0
%     fy =
% 
%              0         0         0         0         0
%        -6.5000   -9.0000    6.0000    6.0000    3.5000
%        -6.5000    3.5000    6.0000    3.5000   -6.5000
%         3.5000    6.0000    6.0000   -9.0000   -6.5000
%              0         0         0         0         0         
%
%     See also GRADIENT, DEL2, AM_LAPLACIAN, AMT, AM_VFC, AM_VFK, AM_GVF,
%     AC_DEFORM, EXAMPLE_VFC, EXAMPLE_PIG. 
% 
%     Reference
%     [1] Bing Li and Scott T. Acton, "Active contour external force using
%     vector field convolution for image segmentation," Image Processing,
%     IEEE Trans. on, vol. 16, pp. 2096-2106, 2007.  
%     [2] Bing Li and Scott T. Acton, "Automatic Active Model
%     Initialization via Poisson Inverse Gradient," Image Processing,
%     IEEE Trans. on, vol. 17, pp. 1406-1420, 2008.   
% 
% (c) Copyright Bing Li 2005 - 2009.

% Revision Log
%   11-30-2005  original 
%   01-20-2006  c version 
%   01-30-2009  minor bug fix

%% inputs and outputs check
if ~ismember(nargin, 1:2) || ndims(A) > 3
    error('Invalid inputs to AM_GRADIENT!')
end

if (ndims(A) == 2 && nargout ~= 2) || (ndims(A) == 3 && nargout ~= 3)   
    error('Invalid outputs to AM_GRADIENT!')
end

%% this part is matlab prototype, please refer to .c and .dll files
if 0,
    if size(A,2)>1,
        varargout{1} = (A(:,[2:end,end-1],:)-A(:,[2,1:end-1],:))/2;
    else
        varargout{1} = 0;
    end
    if size(A,2)>1,
        varargout{2} = (A([2:end,end-1],:,:)-A([2,1:end-1],:,:))/2;
    else
        varargout{2} = 0;
    end
    if ndims(A)==3,
        varargout{3} = (A(:,:,[2:end,end-1])-A(:,:,[2,1:end-1]))/2;
    end
    
    if nargin == 2,
        varargout{1} = varargout{1}/R(1);
        varargout{2} = varargout{2}/R(2);
        varargout{3} = varargout{3}/R(3);        
    end
%% call faster c routine
else
    A = single(A); % using single (32-bits) instead of double (64-bits) precision
    % using c code
    if ndims(A)==2,
        [varargout{1} varargout{2}] = AM_gradient_c(A);
    else
        if nargin == 1,
            [varargout{1} varargout{2} varargout{3}] = AM_gradient_c(A);
        else
            [varargout{1} varargout{2} varargout{3}] = AM_gradient_c(A, single(R));
        end
    end
end