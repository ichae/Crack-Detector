function vert = AM_FFS(Fext,f)
% AM_FFS        Initialize active contours using the FFS algorithm [4].
%     vertices = AM_FFS(Fext, f)
% 
%     Inputs
%     Fext        the external force field. For 2D, d1-by-d2-by-2 matrix, 
%                 the force at (x,y) is [Fext(y,x,1) Fext(y,x,2)]. For 3D,
%                 d1-by-d2-by-d3-by-3 matrix, the force at (x,y,z) is
%                 [Fext(y,x,z,1) Fext(y,x,z,2) Fext(y,x,z,3)].
%     f           binary edge map, where 1's indicate edges. d1-by-d2
%                 matrix for 2D, and d1-by-d2-by-d3 matrix for 3D. 
% 
%     Outputs
%     vertices    vertices position of M initial active models. For 2D,
%                 vertices is a M-by-1 cell vector, where each cell element
%                 is a n-by-2 matrix (n is the number of vertices of this
%                 AC). For 3D, vertices is a M-by-1 structure array, each
%                 element of which contains the faces and vertices and can
%                 be pass directly to the PATCH command.    
% 
%     Note that GRAPHCONNCOMP() in the Bioinformatics Toolbox is required
%     to run this function.
% 
%     Example
%         See EXAMPLE_PIG.
%
%     See also GRAPHCONNCOMP, AMT, AM_PIG, AC_ISOLINE, AM_COD, AC_REMESH, AC_DISPLAY,
%     AM_VFC, AM_VFK, AM_GVF, EXAMPLE_VFC, EXAMPLE_PIG. 
% 
%     Reference
%     [2] Bing Li and Scott T. Acton, "Automatic Active Model
%     Initialization via Poisson Inverse Gradient," Image Processing,
%     IEEE Trans. on, vol. 17, pp. 1406-1420, 2008.   
%     [4] Chunming Li, Jundong Liu, and Martin D. Fox, "Segmentation of
%     external force field for automatic initialization and splitting of
%     snakes," Pattern Recognition, vol.38, pp. 1947–1960, 2005.      
% 
% (c) Copyright Bing Li 2005 - 2009.
 
% Revision Log
%   01-30-2009  original

%% Force field segmentation (FFS) algorithm
% 1. quantize the force field
% 2. convert the quantized force field to a graph using the two conditions 
%    (Q = P+v(P), f(P)>T & f(Q)>T) in the form of sparse adjacency matrix
% 3. use depth-first search algorithm to find weakly connected components
% 4. erode the capture range slightly
% 5. use boundary tracing algorithm to intial snakes

%% inputs check
if ~ismember(nargin, 2) || ~((ndims(Fext) == 3 && size(Fext,3)==2) || (ndims(Fext) == 4 && size(Fext,4)==3)), 
    error('Invalid inputs to AM_FFS!')
end
szF = size(Fext);
if any(size(f) ~= szF(1:end-1))
    error('Invalid inputs!')
end

%% 2D
if ndims(Fext) == 3 
    % 1. quantize the force field
    %      3  2  1
    %       \ | /
    %      4- + -0
    %       / | \
    %      5  6  7
    direct = angle(Fext(:,:,1)-j*Fext(:,:,2))+pi/8;  %[-pi pi]+pi/8
    direct(direct<0) =  direct(direct<0)+2*pi;  %[0 2*pi]
    direct = floor(direct/(pi/4));
    
    % 2. convert the quantized force field to a graph using the two conditions
    %    (Q = P+v(P), f(P)>T & f(Q)>T) in the form of sparse adjacency matrix
    sz = size(Fext);
    sz(end) = [];
    N = prod(sz)+1;
    A = sparse(N,N); % last one represents the infinite point
    for d = 0:7, % convert vector field to directed graph
        [r, c] = find(direct == d);
        idx = sub2ind(sz,r,c);
        switch (d)
            case 0
                B = sparse(idx, mysub2ind(sz,r,c+1), 1,N,N);
            case 1
                B = sparse(idx, mysub2ind(sz,r-1,c+1), 1, N,N);
            case 2
                B = sparse(idx, mysub2ind(sz,r-1,c), 1, N,N);
            case 3
                B = sparse(idx, mysub2ind(sz,r-1,c-1), 1, N,N);
            case 4
                B = sparse(idx, mysub2ind(sz,r,c-1), 1, N,N);
            case 5
                B = sparse(idx, mysub2ind(sz,r+1,c-1), 1, N,N);
            case 6
                B = sparse(idx, mysub2ind(sz,r+1,c), 1, N,N);
            case 7
                B = sparse(idx, mysub2ind(sz,r+1,c+1), 1, N,N);
        end
        A = max(A,B);
    end

    f = f>0;
    for d = 0:3, % convert the edge map to undirected graph
        switch (d)
            case 0
                [r c] = find(f(:,1:end-1) & f(:,2:end));
                idx = sub2ind(sz,r,c);
                B = sparse(idx, sub2ind(sz,r,c+1), 1, N,N);
            case 1
                [r c] = find(f(2:end,1:end-1) & f(1:end-1,2:end));
                idx = sub2ind(sz,r+1,c);
                B = sparse(idx, sub2ind(sz,r,c+1), 1, N,N);
            case 2
                [r c] = find(f(1:end-1,:) & f(2:end,:));
                idx = sub2ind(sz,r,c);
                B = sparse(idx, sub2ind(sz,r+1,c), 1, N,N);
            case 3
                [r c] = find(f(1:end-1,1:end-1) & f(2:end,2:end));
                idx = sub2ind(sz,r,c);
                B = sparse(idx, sub2ind(sz,r+1,c+1), 1, N,N);
        end
        A = max(A,max(B,B'));
    end

    % 3. use depth-first search algorithm to find weakly connected components
    try
        [S, C] = graphconncomp(A,'Weak', true);
    catch
        error('GRAPHCONNCOMP() in the Bioinformatics Toolbox is required to run this function!')
    end
    L = reshape(C(1:end-1), sz);
    ls = unique(C);
    vert = [];

    % 4. erode the capture range slightly
    % 5. use boundary tracing algorithm to intial snakes
    for k = 1:length(ls)
        BW = imerode(L == k,ones(3));
        B = bwboundaries(BW);
        vert = [vert;B];
    end

    for k = 1:length(vert)
        vert{k} = vert{k}(:,[2 1]); 
    end

%% 3D
elseif ndims(Fext) == 4 
    % 1. quantize the force field
    %      3  2  1
    %       \ | /
    %      4- + -0
    %       / | \
    %      5  6  7
    direct = angle(Fext(:,:,:,1)-j*Fext(:,:,:,2))+pi/8;  %[-pi pi]+pi/8
    direct(direct<0) =  direct(direct<0)+2*pi;  %[0 2*pi]
    direct = floor(direct/(pi/4));
    %      0  1
    %      | /
    %      + -2
    %      | \
    %      4  3
    Fmag = sqrt(sum(Fext.^2,4))+eps;
    zenith = acos(Fext(:,:,:,3)./Fmag);
    zenith = floor((zenith+pi/8)/(pi/4));

    % 2. convert the quantized force field to a graph using the two conditions
    %    (Q = P+v(P), f(P)>T & f(Q)>T) in the form of sparse adjacency matrix
    sz = size(Fext);
    sz(end) = [];
    N = prod(sz)+1;
    A = sparse(N,N); % last one represent infinite point
    for z = 0:4
        idx = find(zenith == z);
        [iz, jz, kz] = ind2sub(sz,idx);
        if z == 0,
            B = sparse(idx, mysub2ind(sz,iz,jz,kz+1), 1,N,N);
            A = max(A,B);
        elseif z == 4,
            B = sparse(idx, mysub2ind(sz,iz,jz,kz-1), 1,N,N);
            A = max(A,B);
        else  % z == 1 2 3
            for d = 0:7, % convert vector field to directed graph
                idx = find(zenith == z & direct == d);
                [r, c, kz] = ind2sub(sz,idx);
                switch (d)
                    case 0
                        B = sparse(idx, mysub2ind(sz,r,c+1,kz+2-z), 1,N,N);
                    case 1
                        B = sparse(idx, mysub2ind(sz,r-1,c+1,kz+2-z), 1, N,N);
                    case 2
                        B = sparse(idx, mysub2ind(sz,r-1,c,kz+2-z), 1, N,N);
                    case 3
                        B = sparse(idx, mysub2ind(sz,r-1,c-1,kz+2-z), 1, N,N);
                    case 4
                        B = sparse(idx, mysub2ind(sz,r,c-1,kz+2-z), 1, N,N);
                    case 5
                        B = sparse(idx, mysub2ind(sz,r+1,c-1,kz+2-z), 1, N,N);
                    case 6
                        B = sparse(idx, mysub2ind(sz,r+1,c,kz+2-z), 1, N,N);
                    case 7
                        B = sparse(idx, mysub2ind(sz,r+1,c+1,kz+2-z), 1, N,N);
                end
                A = max(A,B);
            end
        end
    end

    f = f>0;
    for z = -1:1
        for d = 0:3, % convert the edge map to undirected graph
            switch (d)
                case 0
                    offset = [0 1 z];
                case 1
                    offset = [-1 1 z];
                case 2
                    offset = [1 0 z];
                case 3
                    offset = [1 1 z];
            end
            if offset(1) == 1
                i1 = 1:sz(1)-1;  i2 = 2:sz(1);
            elseif offset(1) == -1
                i2 = 1:sz(1)-1;  i1 = 2:sz(1);
            else
                i1 = 1:sz(1);  i2 = i1;
            end
            if offset(2) == 1
                j1 = 1:sz(2)-1;  j2 = 2:sz(2);
            elseif offset(2) == -1
                j2 = 1:sz(2)-1;  j1 = 2:sz(2);
            else
                j1 = 1:sz(2);  j2 = j1;
            end
            if offset(3) == 1
                k1 = 1:sz(3)-1;  k2 = 2:sz(3);
            elseif offset(3) == -1
                k2 = 1:sz(3)-1;  k1 = 2:sz(3);
            else
                k1 = 1:sz(3);  k2 = k1;
            end
            idx = find(f(i1,j1,k1) & f(i2,j2,k2));
            [r c k] = ind2sub(size(f(i1,j1,k1)),idx);
            idx = sub2ind(sz,r-min(0,offset(1)),c-min(0,offset(2)),k-min(0,offset(3)));
            B = sparse(idx, sub2ind(sz,r+max(0,offset(1)),c+max(0,offset(2)),k+max(0,offset(3))), 1, N,N);
            A = max(A,max(B,B'));
        end
    end

    offset = [0 0 1];
    idx = find(f(:,:,1:end-1) & f(:,:,2:end));
    [r c k] = ind2sub(sz-[0 0 1],idx);
    idx = sub2ind(sz-[0 0 1],r-min(0,offset(1)),c-min(0,offset(2)),k-min(0,offset(3)));

    B = sparse(idx, sub2ind(sz,r+max(0,offset(1)),c+max(0,offset(2)),k+max(0,offset(3))), 1, N,N);
    A = max(A,max(B,B'));

    % 3. use depth-first search algorithm to find weakly connected components
    try
        [S, C] = graphconncomp(A,'Weak', true);
    catch
        error('GRAPHCONNCOMP() in the Bioinformatics Toolbox is required to run this function!')
    end
    L = reshape(C(1:end-1), sz);
    ls = unique(C);

    % 4. erode the capture range slightly
    % 5. use boundary tracing algorithm to intial snakes
    ii = 1;
    L = padarray(L,[1 1 1]);  % to form a close surface
    for k = 1:length(ls)
        BW = L == k;
        if length(find(BW)) < N/100
            continue
        end
        BW = imerode(BW,ones([1 1 5]));
        BW = imerode(BW,ones([1 5 1]));
        BW = imerode(BW,ones([5 1 1]));
        fv = isosurface(BW,.5);
        fv.vertices = fv.vertices - 1;  % correct the zero padding
        if ~isempty(fv.vertices)
            if size(fv.vertices,1)>10000
                tmp = reducepatch(fv,10000);%/size(fv(k).vertices,1));
                fv.vertices = tmp.vertices;
                fv.faces = tmp.faces;
            end
            vert(ii) = fv;
            ii = ii+1;
        end
    end
end

%% mysub2ind()
function ind = mysub2ind(sz, i, j, k)
N = prod(sz);
ind = (N+1)*ones(size(i));
if nargin == 3
    inrange = 1<=i & i<=sz(1) & 1<=j & j<=sz(2);
    tmp = (j-1)*sz(1)+i;
else
    inrange = 1<=i & i<=sz(1) & 1<=j & j<=sz(2) & 1<=k & k<=sz(3);
    tmp = ((k-1)*sz(2)+j-1)*sz(1)+i;
end
ind(inrange) = tmp(inrange);