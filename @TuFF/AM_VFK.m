function K = AM_VFK(D, r, type, a, R)
% AM_VFK        Compute the vector field kernel (VFK).
%     K = AM_VFK(D, r, 'power', alpha)
%     K = AM_VFK(..., R)
%   
%     Inputs
%     D           output dimension, 2 for 2D VFK, 3 for 3D VFK.
%     r           the VFK radius.
%     gamma       decrease parameter of equation (15) in [1].
%     R           a 3-element vector, [Rx Ry Rz] define the spacing between 
%                 voxels for x, y and z dimension, default value is [1 1 1].
%               
%     Outputs
%     K           the VFC external force kernel, d=2*r+1. For 2D,
%                 d-by-d-by-2 matrix, the force at (x,y) is [K(y,x,1)
%                 K(y,x,2)]. For 3D, d-by-d-by-d-by-3 matrix, the force at
%                 (x,y,z) is [K(y,x,z,1) K(y,x,z,2) K(y,x,z,3)].  
% 
%     Note that for memory saving reason, the output data class is single / float (32-bit).
%     
%     Example
%         K = AM_VFK(2, 4, 'power', 2);
%         AC_QUIVER(K);
%
%     See also AMT, AM_VFC, AM_PIG, AM_ISOLINE, AM_DEFORM, AC_INITIAL, AC_REMESH,
%     AC_DISPLAY, AM_GVF, EXAMPLE_VFC, EXAMPLE_PIG. 
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
%   02-11-2005  original
%   01-30-2009  minor bug fix

%% inputs check
if ~ismember(nargin, 4:5) || r<1 || ~ismember(D, 2:3)
    error('Invalid inputs to AM_VFK!')
elseif nargin == 4
    R = ones(1,D);
end

% if ~strcmp(type, 'power')
%     error('Only support power magnitude function right now!')
% end
%%
r0 = floor(r./R);
if D == 2
    [x y] = meshgrid(R(1)*(r0(1):-1:-r0(1)),R(2)*(r0(2):-1:-r0(2)));
    dist = sqrt(x.*x+y.*y);
    SetZero = (dist>r);
    x(SetZero) = 0;
    y(SetZero) = 0;

    if strcmp(type, 'power')
        dist(~SetZero) = dist(~SetZero).^(a+1);
        dist((end+1)/2) = 1e-8; % set the distance at the center a small value to prevent division by zero
        x(~SetZero) = x(~SetZero)./dist(~SetZero);
        y(~SetZero) = y(~SetZero)./dist(~SetZero);
    elseif strcmp(type, 'gaussian')
        dist(~SetZero) = exp(dist(~SetZero).^2/(a*a));
        dist((end+1)/2) = 1e-8; % set the distance at the center a small value to prevent division by zero
        x(~SetZero) = (1/(2*a*a))*x(~SetZero)./dist(~SetZero);
        y(~SetZero) = (1/(2*a*a))*y(~SetZero)./dist(~SetZero);
    end
    
    K = cat(3,x,y);
elseif D == 3
    [x y z] = meshgrid(R(1)*(r0(1):-1:-r0(1)),R(2)*(r0(2):-1:-r0(2)),R(3)*(r0(3):-1:-r0(3)));
    dist = sqrt(x.*x+y.*y+z.*z);
    SetZero = (dist>r);
    x(SetZero) = 0;
    y(SetZero) = 0;
    z(SetZero) = 0;

    if strcmp(type, 'power')
        dist(~SetZero) = dist(~SetZero).^(a+1);
        dist((end+1)/2) = 1e-8; % set the distance at the center a small value to prevent division by zero
        x(~SetZero) = x(~SetZero)./dist(~SetZero);
        y(~SetZero) = y(~SetZero)./dist(~SetZero);
        z(~SetZero) = z(~SetZero)./dist(~SetZero);
    end
    K = single(cat(4,x,y,z));
else
    error('D must be 2 or 3!');   
end
