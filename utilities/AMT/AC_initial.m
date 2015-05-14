function vertex = AC_initial(res, type, p)
% AC_INITIAL     Initialize an active contour (AC), also known as snake.
%     vertex = AC_INITIAL(res)
%     Select points from the figure using mouse to form a close contour with 
%     resolution (distance between vertices) res. See getline() for instructions.
% 
%     vertex = AC_INITIAL(res, 'open')
%     Select points from the figure using mouse to form an open contour.
% 
%     vertex = AC_INITIAL(res, 'circle', [cx cy r])
%     Initialize a circle centered at [cx cy] with radius r. 
% 
%     vertex = AC_INITIAL(res, 'rect', [xmin xmax ymin ymax])
%     Initialize a rectangle whose top-left corner at [xmin ymin]
%     and bottom-right corner at [xmax ymax]
% 
%     Example
%         I = imread('im_U.bmp');
%         vert = AC_initial(.5, 'rect', [10 50 20 40]);
%         imshow(I);
%         AC_display(vert);
%
%     See also GETLINE, AMT, AC_REMESH, AC_DISPLAY, AM_VFC, AM_VFK, AM_PIG, 
%     EXAMPLE_VFC, EXAMPLE_PIG. 
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
%   01-30-2009  minor bug fix

%% inputs check
if ~ismember(nargin, 1:3)
    error('Invalid inputs to AC_INITIAL!')
end
if nargin == 2 && ~strcmp(type, 'open')
    error('Invalid inputs to AC_INITIAL!')
end

%% initialize
if nargin <3,    
    [x y] = getline;
else
    if strcmp(type, 'circle')
        if numel(p) == 3
            t = (pi/10:pi/10:pi*2)';
            x = p(1) + p(3)*cos(t);
            y = p(2) + p(3)*sin(t);
        else
            error('Invalid inputs to AC_INITIAL!')
        end
    elseif strcmp(type, 'rect')
        if numel(p) == 4
            x = p([1 2 2 1]).';
            y = p([3 3 4 4]).';
        else
            error('Invalid inputs to AC_INITIAL!')
        end
    else
        error('Invalid type!')
    end
end
%% remesh to desired resolution
if nargin == 2, % open contour
    vertex = AC_remesh([x y], res, 'adaptive', 'open');
else
    vertex = AC_remesh([x y], res);
end