function [vert COD] = AM_CoD(Fext,R)
% AM_COD        Initialize active contours using the CoD algorithm [3].
%     [vertices CoD] = AM_COD(Fext, R)
% 
%     Inputs
%     Fext        the external force field. For 2D, d1-by-d2-by-2 matrix, 
%                 the force at (x,y) is [Fext(y,x,1) Fext(y,x,2)]. For 3D,
%                 d1-by-d2-by-d3-by-3 matrix, the force at (x,y,z) is
%                 [Fext(y,x,z,1) Fext(y,x,z,2) Fext(y,x,z,3)].
%     R           if the distance between centers of divergence (CoD) is
%                 less than R, they are considered as one source. One small
%                 circle with radius R is place on each source. 
%     Outputs
%     vertices    vertices position of M initial active models. For 2D,
%                 vertices is a M-by-1 cell vector, where each cell element
%                 is a n-by-2 matrix (n is the number of vertices of this
%                 AC). For 3D, vertices is a M-by-1 structure array, each
%                 element of which contains the faces and vertices and can
%                 be pass directly to the PATCH command.    
%     CoD         binary d1-by-d2 matrix, where 1's indicate the CoD locations.
% 
%     Example
%         See EXAMPLE_PIG.
%
%     See also AMT, AM_PIG, AC_ISOLINE, AM_FFS, AC_REMESH, AC_DISPLAY,
%     AM_VFC, AM_VFK, AM_GVF, EXAMPLE_VFC, EXAMPLE_PIG. 
% 
%     Reference
%     [2] Bing Li and Scott T. Acton, "Automatic Active Model
%     Initialization via Poisson Inverse Gradient," Image Processing,
%     IEEE Trans. on, vol. 17, pp. 1406-1420, 2008.   
%     [3] Xingfei Ge and Jie Tian, "An automatic active contour model for
%     multiple objects," in Pattern Recognition, Proceedings of the
%     International Conference on, 2002.    
% 
% (c) Copyright Bing Li 2005 - 2009.
 
% Revision Log
%   01-30-2009  original

%% inputs check
if ~ismember(nargin, 2) || ~((ndims(Fext) == 3 && size(Fext,3)==2) || (ndims(Fext) == 4 && size(Fext,4)==3)) || numel(R) ~= 1, 
    error('Invalid inputs to AM_COD!')
end

if ndims(Fext) == 3, % 2D 
    x = Fext(:,:,1);
    x1 = x(:,[2:end end]);    
    Psx = x<=x1 & abs(sign(x)+sign(x1))<=1;
    
    y = Fext(:,:,2);
    y1 = y([2:end end],:);
    Psy = y<=y1 & abs(sign(y)+sign(y1))<=1;
    
    COD = Psx & Psy;
    
    SE = strel('disk',floor(R/2));
    mask = imdilate(COD,SE);    % merge CoD that are close to each other
    STATS = regionprops(bwlabel(mask),'Centroid');  % find where to put the circles
    for k = 1:length(STATS)
        vert{k} = AC_initial(1,'circle', [STATS(k).Centroid R]);
    end    
elseif ndims(Fext) == 4, % 3D 
    error('AM_CoD does not support 3D vector field right now!');
    
    COD = Fext(:,:,:,1) <= Fext(:,[2:end end],:,1) &  ...                  % x<=x1
        abs(sign(Fext(:,:,:,1))+sign(Fext(:,[2:end end],:,1)))<=1  &  ...  % abs()<=1      
        Fext(:,:,:,2) <= Fext([2:end end],:,:,2) &        ...              % y<=y1
        abs(sign(Fext(:,:,:,2))+sign(Fext([2:end end],:,:,2)))<=1  & ...   % abs()<=1      
        Fext(:,:,:,3) <= Fext(:,:,[2:end end],3) &        ...              % z<=z1
        abs(sign(Fext(:,:,:,3))+sign(Fext(:,:,[2:end end],3)))<=1;  ...    % abs()<=1      
    
    for k = 1:3
        sz = ones(1,3);
        sz(k) = floor(R/2);
        SE = ones(sz);
        mask = imdilate(COD,SE);    % merge COD that are close to each other
    end
    L = bwlabeln(mask);
    for k = max(L(:)):-1:1
        idx = find(L == k);
        [y, x, z] = ind2sub(size(L),idx);
        vert(k) = AS_initial(1,'sphere', [mean(x) mean(y) mean(z) R]);   % AS_initial is under development
    end    
else
    error('wrong input, only support 2D and 3D now!')
end
     
