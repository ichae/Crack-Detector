function E = AM_PIG(Fext, omega, f)
% AM_PIG        Estimate the energy via Poisson inverse gradient (PIG) [2].
%     E = AM_PIG(Fext, omega, f)
% 
%     Inputs
%     Fext        the external force field. For 2D, d1-by-d2-by-2 matrix, 
%                 the force at (x,y) is [Fext(y,x,1) Fext(y,x,2)]. For 3D,
%                 d1-by-d2-by-d3-by-3 matrix, the force at (x,y,z) is
%                 [Fext(y,x,z,1) Fext(y,x,z,2) Fext(y,x,z,3)].
%     omega       binary mask, where 1's indicate pixels to be
%                 interpolated. d1-by-d2 matrix for 2D, and d1-by-d2-by-d3
%                 matrix for 3D.   
%     f           edge map, the pixel values at the outter boundary of
%                 omega is used in the Dirichlet boundary condition.
%                 d1-by-d2 matrix for 2D, and d1-by-d2-by-d3 matrix for 3D.   
% 
%     Outputs
%     E           the external estimate energy calculated from the external
%                 force field via PIG. Note that E(~omega) = -f(~omega).
% 
%     Example
%         See EXAMPLE_PIG.
%
%     See also AMT, AC_ISOLINE, AM_COD, AM_FFS, AC_REMESH, AC_DISPLAY,
%     AM_VFC, AM_VFK, AM_GVF, EXAMPLE_VFC, EXAMPLE_PIG. 
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
 
% Revision Log:
%   09-17-2006	original
%   01-30-2009  minor bug fix

%% inputs check
if nargin ~= 3
    error('Invalid input!')
end

if (ndims(Fext) == 4 && size(Fext,4)==3) 
    Dim = 3;    % 3D
elseif (ndims(Fext) == 3 && size(Fext,3)==2)
    Dim = 2;    % 2D
else
    error('Invalid input: incorrect Fext size!')
end

if ~isa(omega,'logical')
    error('Invalid input: omega must be logical!')
end

szFext = size(Fext);
szFext = szFext(1:Dim);  % a-by-b or a-by-b-by-c
if any(size(f)~=size(omega)) || any(size(f)~=szFext) || any(size(omega)~=szFext)
    error('Invalid input: the input size must match!')
end

%% Equation (18)
E = -f;     

%% Solve Equation (19) in matrix form - A*E = b*f + div(Fext) = B
pixnum = prod(szFext);	% # of total pixels, a*b or a*b*c
idx = find(omega);      % find the pixels to be interpolated
N = length(idx);        % # of pixels to be interpolated
if N==0     % nothing to interpolate
    return
end
[ti tj tk] = ind2sub(szFext,idx);
t = cat(2, ti, tj, tk);

%% generate A, representing the relationship between pixels in omega
A = sparse(1:N, 1:N, -Dim*2);
for i = 1:Dim,
    for j = [-1 1]
        tmp = t;
        tmp(:,i) = tmp(:,i)+j;
        [tf, neighbor] = ismember(tmp, t, 'rows'); % check if the neighbor is in omega
        A = A + sparse(find(tf), neighbor(tf) , 1, N, N);
    end
end

% take care of image boundaries (first or last row/colomn)
for i = 1:Dim,
    tind = find(t(:,i) == 1 | t(:,i) == szFext(i));
    A = A + sparse(tind, tind, 1, N, N);
end

%% generate b, representing the boundary condition b*f(Bidx)
Bidx = find(bwperim(~omega));   % outter boundary pixels
NB = length(Bidx);              % # of boundary pixels
b = sparse(N, NB);              % initialize b
[ti tj tk] = ind2sub(szFext,Bidx);
Bt = cat(2, ti, tj, tk);
clear ti tj tk tind

for i=1:Dim,
    for j = [-1 1]
        tmp = t;
        tmp(:,i) = tmp(:,i)+j;
        [tf, neighbor] = ismember(tmp, Bt, 'rows');
        b = b + sparse(find(tf), neighbor(tf), 1, N, NB);
    end
end

%% generate B = b*f + div(Fext)
B = b*double(f(Bidx));
clear b Bidx neighbor Bt tf

for i = 1:Dim
    for j = [-1 1]
        tmp = t;
        tmp(:,i) = tmp(:,i)+j;
        i1 = tmp(:,i) == 0 | tmp(:,i)==szFext(i)+1;
        tmp = trunc(szFext, tmp(:,1:Dim));

        if Dim == 2
            ind = sub2ind([szFext 2], tmp(:,1), tmp(:,2), (3-i)*ones(size(tmp,1),1));
        else
            ind = sub2ind([szFext 3], tmp(:,1), tmp(:,2), tmp(:,3), (mod(2-i,3)+1)*ones(size(tmp,1),1));
        end
        B = B - 0.5*j*Fext(ind).*(-1).^i1;  % i1 takes care of div() on image boundaries
    end
end

%% solve A*E = B to interpolate E within omega
clear Fext t tmp ind i1 % clear to save memory for the next step
E(idx) = A\double(B);

%% function sub = trunc(sz, sub);
function sub = trunc(sz, sub);
% make sure the index or subscript is within the valid range, truncate
%  
%  sub = trunc(sz, sub);
%       sz  -  size of the range, [a b] or [a b c]
%       sub -  subscripts, n-by-2 or 3 array, organized as follow
%               [y1 x1 z1]
%               [y2 x2 z2]
%               [. . . . ]  
%               where y and x are row and column index, respectively

if nargin ~=2 || length(sz) ~= size(sub,2) || size(sub,2)>3
    error('Invalid input!')
end

for i = 1:length(sz)
    sub(:,i) = min(max(1,sub(:,i)),sz(i));
end