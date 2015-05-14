function B = AM_laplacian(A)
% AM_LAPLACIAN    2D/3D discrete Laplacian, a faster implementation than build-in DEL2().
%     L = AM_laplacian(U) where U is a 2-D or 3-D matrix
% 
%     The output L has the same size as the input U.
% 
%     Note that the results at the boundaries are different from the ones
%     generated from DEL2().    
% 
%     Example
%     lap = AM_laplacian(magic(5))
%    
%     lap =
%         6.5000  -17.0000   10.5000    3.0000   -3.0000
%       -15.2500   10.0000    1.2500   -1.2500    0.2500
%         7.2500    2.5000         0   -2.5000   -7.2500
%        -0.2500    1.2500   -1.2500  -10.0000   15.2500
%         3.0000   -3.0000  -10.5000   17.0000   -6.5000
%
%     See also DEL2, GRADIENT, AM_GRADIENT, AMT, AM_VFC, AM_VFK, AM_GVF,
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
if ~ismember(nargin, 1) || ndims(A) > 3
    error('Invalid inputs to AM_LAPLACIAN!')
end

%% this part is matlab prototype, please refer to .c and .dll files
if 0,
    if ndims(A)==3,
        B = (A([2,1:end-1],:,:) + A([2:end,end-1],:,:) + A(:,[2,1:end-1],:) + ...
            A(:,[2:end,end-1],:) + A(:,:,[2,1:end-1]) + A(:,:,[2:end,end-1]))/6-A;
    else
        B = (A([2,1:end-1],:) + A([2:end,end-1],:) + A(:,[2,1:end-1]) + A(:,[2:end,end-1]))/4 - A;
    end
%% call faster c routine
else
    B = AM_laplacian_c(double(A));
end