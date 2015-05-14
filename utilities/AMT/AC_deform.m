function vertex = AC_deform(vertex,alpha,beta,tau,Fext,ITER,type)
% AC_DEFORM     Deform an active contour (AC), also known as snake.
%     vertex1 = AC_DEFORM(vertex0,alpha,beta,tau,Fext,iter)
%     vertex1 = AC_DEFORM(vertex0,alpha,beta,tau,Fext,iter,type)
%   
%     Inputs
%     vertex0     position of the vertices, n-by-2 matrix, each row of 
%                 which is [x y]. n is the number of vertices.
%     alpha       AC elasticity (1st order) parameter ranges from 0 to 1.
%     beta        AC rigidity (2nd order) parameter ranges from 0 to 1.
%     tau         time step of each iteration.
%     Fext        the external force field,d1-by-d2-by-2 matrix, 
%                 the force at (x,y) is [Fext(y,x,1) Fext(y,x,2)].
%     iter        number of iterations, usually ranges from 1 to 5.
%     type        'close' - close contour (default), the last vertex and
%                           first vertex are connected 
%                 'open'  - open contour,  the last vertex and first vertex
%                           are not connected 
%               
%     Outputs
%     vertex1     position of the vertices after deformation, n-by-2 matrix
%     
%     Note that if the vertices are outside the valid range, i.e., y>d1 ||
%     y<1 || x>d2 || x<1, they will be pulled inside the valid range. 
% 
%     Example
%         See EXAMPLE_VFC, EXAMPLE_PIG.
%
%     See also AMT, AM_VFC, AM_VFK, AM_PIG, AC_INITIAL, AC_REMESH,
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
%   01-30-2006  external force interpolation outside the image 
%   02-18-2006  add open contour codes
%   01-30-2009  minor bug fix

%% inputs check
if ~ismember(nargin, 6:7) || ndims(Fext) ~= 3 || size(Fext,3) ~= 2,
    error('Invalid inputs to AC_DEFORM!')    
end
if nargin == 6
    type = 'close';
end

N = size(vertex,1);
if size(vertex,2) ~= 2
    error('Invalid vertex matrix!')
end

if N < 3
    return
end

%% compute T = (I + tao*A) of equation (9) in reference [1]
Lap = sparse(1:N, 1:N, -2) + sparse(1:N, [N 1:N-1], 1) + sparse(1:N, [2:N 1], 1);

if strcmp(type,'open'), % offset tau for boundary vertices    
    tau = sparse(1:N, 1:N, tau);
    tau(1) = 0;
    tau(end) = 0;
    offset = sparse(1:N, 1:N, 1);
    offset(1,1)=0;      offset(N,N) = 0;
    offset(1,2)=1;      offset(N,N-1) = 1;
    Lap = offset*Lap;
end

T =sparse(1:N,1:N,1)+ tau*(beta*Lap*Lap-alpha*Lap);

%% Another way to compute T for close AC
% a = beta;
% b = -alpha - 4*beta;
% c = 2*alpha + 6*beta;
% 
% T = sparse(1:N,1:N,1) + tau*(sparse(1:N,1:N,c) + sparse(1:N,[N,1:N-1],b) + sparse(1:N,[2:N,1],b)...
%     + sparse(1:N,[N-1,N,1:N-2],a) + sparse(1:N,[3:N,1,2],a));

%% Deform
center = size(Fext)/2;
center = center([2,1]);
for i=1:ITER,
    IdxOut = find(vertex(:,1)<1 | vertex(:,1)>size(Fext,2) | vertex(:,2)<1 | vertex(:,2)>size(Fext,1));
    IdxIn = setdiff(1:N,IdxOut);
    for i=1:2,
        % interpolate the external force for vertices within the range
        F(IdxIn,i)  = interp2(Fext(:,:,i),vertex(IdxIn,1),vertex(IdxIn,2),'*linear');
        F(IdxOut,i) = center(i)-vertex(IdxOut,i);  % pointing to the image center
    end
    if ~isempty(IdxOut)
        % normalize the forces outside the image
        Fmag = sqrt(sum(F(IdxOut,:).^2,2));
        F(IdxOut,1) = F(IdxOut,1)./Fmag;
        F(IdxOut,2) = F(IdxOut,2)./Fmag;
    end       

    vertex = T \(vertex+tau*double(F));   % equation (9) in reference [1]
end