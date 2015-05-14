function vertex = AC_remesh(vertex,res,method, type)
% AC_REMESH     Remesh an active contour (AC) to ensure the resolution.
%     vertex1 = AC_REMESH(vertex0)
%     vertex1 = AC_REMESH(vertex0,res)
%     vertex1 = AC_REMESH(vertex0,res,method)
%     vertex1 = AC_REMESH(vertex0,res,method,type)
% 
%     Inputs
%     vertex0     position of the vertices, n-by-2 matrix, each row of 
%                 which is [x y]. n is the number of vertices.
%     res         desired resolution (distance between vertices) in pixels,
%                 default value is 1.
%     method      'adaptive' - (default) the distances between vertices
%                           will be larger than res/2 and less than res*2,
%                           the shape will be the same
%                 'equal' - the distances between vertices will be the same, 
%                           the shape may be slightly different after remeshing
%     type        'close' - close contour (default), the last vertex and
%                           first vertex are connected 
%                 'open'  - open contour,  the last vertex and first vertex
%                           are not connected 
%     Outputs
%     vertex1     position of the vertices after remeshing, m-by-2 matrix
%     
%     Example
%         See EXAMPLE_VFC, EXAMPLE_PIG.
%
%     See also AMT, AC_DEFORM, AM_VFC, AM_VFK, AM_PIG, AC_INITIAL, 
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
%   11-30-2005  original 
%   02-18-2006  support open type
%   01-30-2009  minor bug fix

%% inputs check
if ~ismember(nargin, 1:4) || numel(res) ~= 1,
    error('Invalid inputs to AC_REMESH!')    
end
if size(vertex,2) ~= 2 
    error('Invalid vertex matrix!')
end

if nargin<2,
    res = 1;
end    
if nargin<3,
    method = 'adaptive';    
end
if nargin<4,
    type = 'close';    
end
%% delete the vertices close to each other (distance less than res/2)
while 1,
    N = size(vertex,1);
    if N<5,
        return
    end
    p = vertex(:,1)+sqrt(-1)*vertex(:,2);    
    dist = abs(p(:,ones(1,N))-p(:,ones(1,N)).');   % calculate distances between all the vertices
    [i,j] = find(dist<res/2);
    i = i(i~=j);    
    idx = medfilt2(ismember(1:N,i),[1 3]);
    i = union(find(idx(2:2:end))*2,find(imopen(idx,ones(1,5))));
   
    if isempty(i),
        break           % no more close vertices left
    end
    
    try 
        vertex(i,:) = [];
    catch
        break
    end
end

%% remesh
if strcmp(method,'adaptive'),      
    if strcmp(type,'open'),  
        d = sqrt(sum((vertex(2:end,:) - vertex(1:end-1,:)).^2,2));
    else
        d = sqrt(sum((vertex([2:end 1],:) - vertex).^2,2));
    end
    while max(d)>res*2, % insert vertices between faraway vertices
        id = find(d>res*2);
        N = size(vertex,1);
        ind = [1:N,id'+.5];
        vertex = [vertex; (vertex(id,:)+vertex(mod(id,N)+1,:))/2];
        [st,ix] = sort(ind);
        vertex = vertex(ix,:);
        if strcmp(type,'open'),
            d = sqrt(sum((vertex(2:end,:) - vertex(1:end-1,:)).^2,2));
        else
            d = sqrt(sum((vertex([2:end 1],:) - vertex).^2,2));
        end   
    end
elseif strcmp(method,'equal'),
    dif = diff(vertex);    
    d = sqrt(sum(dif.^2,2));  % compute the distance from previous vertex for point 2:N+1
    D = sum(d);
    M = ceil(D/res);
    di = linspace(0,D,M+1);
    dcum = [0;cumsum(d)];

    for i=1:2,
        vert(:,i) = interp1(dcum,vertex(:,i),di(1:end-1))';
    end    
    vertex = vert;
else
    error('Invalid method!');
end