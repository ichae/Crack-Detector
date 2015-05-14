function H = AC_quiver(Fext, I)
% AC_QUIVER    Display the 2D external force field in arrows.
%   AC_QUIVER(Fext)
%   AC_QUIVER(Fext, I)
%   h = AC_QUIVER(...)
%
%     Inputs
%     Fext        the external force field,d1-by-d2-by-2 matrix, 
%                 the force at (x,y) is [Fext(y,x,1) Fext(y,x,2)].
%     I           display the vector field on top of the image I, d1-by-d2 matrix
% 
%     Outputs
%     h           handles to the arrows
% 
%     Example
%         I = imread('im_U.bmp');
%         K = AM_VFK(2, 32, 'power',1.8);
%         Fext = AM_VFC(~I, K, 1);
%         AC_quiver(Fext,I);
%
%     See also QUIVER, AMT, AM_VFC, AM_VFK, AM_GVF, AC_INITIAL, AC_REMESH,
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
%   02-11-2006  original
%   01-14-2009  support AC_quiver(Fext, I) 
%   01-30-2009  minor bug fix

%% inputs check
H = NaN;
if ~ismember(nargin, 1:2) || ndims(Fext) ~= 3 || size(Fext,3) ~= 2,
    error('Invalid inputs to AC_QUIVER!')
end

if nargin == 2
    if size(Fext,1) ~= size(I,1) || size(Fext,2) ~= size(I,2) 
        error('The size of Fext does not match the size of I!')
    end
end

%% quiver
if nargin == 1
    H = quiver(Fext(:,:,1), Fext(:,:,2), .5);
else
    imshow(I);
    colormap gray;
    hold on
    H = quiver(Fext(:,:,1), Fext(:,:,2), .5);
    hold off  
end  
axis ij
axis equal
axis tight
axis off
