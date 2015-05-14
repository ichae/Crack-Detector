function [vert isoValue] = AC_isoLine(E, lambda, eta)  
% AC_ISOLINE    Select the optimal close isoline(s).
%     [vertices isoValue] = AC_ISOLINE(E, lambda)
%     [vertices isoValue] = AC_ISOLINE(E, lambda, eta)
%   
%     Inputs
%     E           AC external energy, d1-by-d2 matrix.
%     lambda      isovalues which isolines are computed at. 
%     eta         the number of isolines to be selected, 1 by default.
%               
%     Outputs
%     vertices    vertices position of eta selected optimal isolines,
%                 eta-by-1 cell vector, where each cell element is a n-by-2                 
%                 matrix (n is the number of vertices of this isoline).
%                 Special case - if eta is 1, then vert is a n-by-2 matrix.
%     isoValue    the corresponding isoValue of the optimal isolines.    
%     
%     Note that this function can be modified for open isoline
%     y<1 || x>d2 || x<1, they will be pulled inside the valid range. 
% 
%     Example
%         See EXAMPLE_PIG.
%
%     See also CONTOURC, AMT, AM_VFC, AM_VFK, AM_PIG, AC_INITIAL, AC_REMESH,
%     AC_DISPLAY, AC_DEFORM, EXAMPLE_VFC, EXAMPLE_PIG. 
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
%   01-30-2009  original

%% inputs check
if ~ismember(nargin, 2:3) || ndims(E) ~= 2 || numel(eta) ~= 1
    error('Invalid inputs to AC_ISOLINE!')    
end

if nargin == 2
    eta = 1;
end
%% select optimal isolines in terms of minimum external energy
vert = [];
isoValue = lambda;
MinExtEng = Inf;
Emax = max(E(:));
for isovalue  = lambda,
    [v leng] = LongestCloseLine(E, isovalue, eta);

    epsilon = 0;
    % equation (25) in reference [2], assume alpha = beta = 0
    for i = 1:length(v)
        epsilon = epsilon + (isovalue-Emax)*leng(i);    
    end
    % find minimum
    if epsilon < MinExtEng,
        MinExtEng = epsilon;
        vert = v;
        isoValue = isovalue;
    end
end
%% special case
if length(v) == 1
    v = v{1};
end


%% %%%%%%%%%%%%%%%%%%
%% LongestCloseLine()
function [v leng] = LongestCloseLine(f, value, k)
% LONGESTCLOSELINE
%     v = LongestCloseLine(f, iv)
%     v = LongestCloseLine(f, iv, k)
%     [v l] = LongestCloseLine(..)
%     Pick the k longest lines given the function f and the isovalue iv,
%     return k-by-1 cell matrix v, and corresponding k-by-1 length vector l.
% 
% (c) Copyright Bing Li 2005 - 2009.

%% inputs check
if ~ismember(nargin, 2:3)
    error('Invalid inputs to LongestCloseLine!')
end
if nargin == 2
    k = 1;
end

%% compute isolines
C = contourc(f,value*[1 1]);
if isempty(C)
    v = []; l = [];  return;
end

%% convert C to v
i = 1;
while ~isempty(C)
    if size(C,2) < C(2,1)+1
        error('incorrect input!')
    end
    
    leng(i) = C(2,1);
    v{i} = C(:,2:C(2,1)+1).';
    C(:,1:C(2,1)+1) = [];
    i = i+1;    
end

%% find close contours
for i = 1:length(v)
    flag(i) = 2>sqrt(sum((v{i}(1,:) - v{i}(end,:)).^2));    % close lines only
end
v = v(flag);
leng = leng(flag);

%% Pick the k longest lines
if length(leng)>k
    [sorted_leng, idx] = sort(leng,2,'descend');
    v = v(idx(1:k));
    leng = leng(idx(1:k));
end
