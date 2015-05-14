function h = AC_display(vertex,flag,LineSpec)
% AC_DISPLAY    Display the contour on the current axis.
%     AC_DISPLAY(vertex)
%     AC_DISPLAY(vertex,type)
%     AC_DISPLAY(vertex,type,LineSpec)
%     h = AC_DISPLAY(...)
%
%     Inputs
%     vertex      position of the vertices, n-by-2 matrix, each row of 
%                 which is [x y]. n is the number of vertices.
%     type        'close' - close contour (default), the last vertex and
%                           first vertex are connected 
%                 'open'  - open contour,  the last vertex and first vertex
%                           are not connected 
%     LineSpec    a line specification that determines line type, marker symbol, 
%                 and color of the plotted lines. Type 'doc LineSpec' for more details.
% 
%     Outputs
%     h           the handle to the contour plotted
% 
%     Example
%         I = imread('im_U.bmp');
%         vert = AC_initial(.5, 'circle', [32 32 24]);
%         imshow(I);
%         h = AC_display(vert,'close','r--');
%         set(h,'LineWidth',3)
%
%     See also PLOT, AMT, AM_VFC, AM_VFK, AM_PIG, AC_INITIAL, AC_REMESH,
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
%   01-30-2009  minor bug fix

%% inputs check
h= NaN;
if ~ismember(nargin, 1:3)
    error('Invalid inputs to AC_DISPLAY!')    
end
if nargin<2,
    flag = 'close';
end    

if size(vertex,1) == 0
    return
end
if size(vertex,2) ~= 2 
    error('Invalid vertex matrix!')
end

if strcmp(flag,'close'),
    vertex = vertex([1:end,1],:);
end

%% plot
hold on
if nargin == 3
   h = plot(vertex(:,1),vertex(:,2),LineSpec);
else
   h = plot(vertex(:,1),vertex(:,2));
end
hold off