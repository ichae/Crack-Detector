function [u,shadow] = shadowRemove2(param)
% Shadow and removal using variational algorithm

f       = param.Img;
type    = param.type;
s       = param.sigma;
n_poly  = param.n_poly;
n_iter  = param.n_iter;
lambda  = param.lambda;
error   = param.error;
count_lim   = param.cnvg_cnt;
interval    = param.intrvl;

% Compute the shadow function
if strcmp(type,'gaussian')
    sz = 6*s;
    h = fspecial('gaussian',sz,s);
    shadow = imfilter(f,h,'same','replicate');
elseif strcmp(type,'none')
    shadow = 1;
elseif strcmp(type,'legendre')
    shadow = findShadow(f,n_poly);
elseif strcmp(type,'minimax')
    shadow = findShadowMinimax(f,300,count_lim,error,s);
end

 
%% Denoising with shadow removal 

J = zeros(n_iter,1);
data_coeff = 1;
u = f;
figure('name','Shadow Removal+Denoising');
subplot(1,3,1);imshow(f);
subplot(1,3,2);imshow(u);
subplot(1,3,3);imshow(shadow);
count = 0;

for its = 1 : n_iter
    u = NeumanBoundCond(u);
    [ux,uy] = gradient(u);
    g = sqrt(ux.^2+uy.^2+eps);
    [uxx,~] = gradient(ux./g);
    [~,uyy] = gradient(uy./g);
    
    term1 = (f-u.*shadow).*shadow;
    term2 = uxx+uyy;
    
    term1 = term1./(max(abs(term1(:))));
    term2 = term2./(max(abs(term2(:))));
    
    chng = data_coeff*term1+ lambda*term2;
    dt = 0.01/max(abs(chng(:)));
    u = u + dt*chng;
    if mod(its,interval) == 0
        subplot(1,3,2);imshow(u/max(u(:))); drawnow;
        subplot(1,3,3);imshow(shadow); drawnow;
        disp(its);
    end
    
    % Objective function
    data_trm = 0.5*sum((f(:)-shadow(:).*u(:)).^2);
    smooth_trm = 0.5*sum(g(:));
    
    J(its) = data_coeff*data_trm + lambda*smooth_trm;
    if its > 20
        count = convergence(J,its,count,error);
    end
    if  count > count_lim
        break;
    end
    
end
u = u/max(u(:));
fprintf('Converged at %d\n',its);
% figure(100);
% newJ = J(1:its);
% plot(1:its,newJ,'-');drawnow; title('convergence');
% shadow = shadow/max(shadow(:));
end

% =========================================================================
% =========================================================================

function s = findShadowMinimax(f,n_iter,count_lim,error,s)
% Estimate the shadow by variational minimax

sz = 6*s;
h1 = fspecial('gaussian',sz,s);
s = imfilter(f,h1,'same','replicate');
% s = f;

p = 2.4;
[gx,gy] = gradient(medfilt2(f));
g = (sqrt(gx.^2+gy.^2)).^p; 
g = g/max(g(:));
% g = 1;
figure('name','Shadow Function');imshow(s);
J = zeros(n_iter,1);
count = 0;
for ii = 1 : n_iter
    s = NeumanBoundCond(s);
    [sx,sy] = gradient(s);
    [sxx,~] = gradient(sx); [~,syy] = gradient(sy);
    
    E1 = (f-s).^2.*g;   
    E1 = sum(E1(:));
    E2 = sx.^2+sy.^2;   
    E2 = sum(E2(:));
    alpha = E2/sqrt(E1.^2+E2.^2);
    term1 = (f-s).*g;
    term2 = (sxx+syy);

    term1 = term1/(max(abs(term1(:))));
    term2 = term2/(max(abs(term2(:))));
    
    chng = sqrt(1-alpha^2)*term1 + alpha*term2;
    dt = 0.01/max(abs(chng(:)));
    s = s + dt*chng;
    disp(alpha);
    imshow(s);drawnow;
    
    % Check convergence
    J(ii) = sqrt(1-alpha^2)*E1 + alpha*E2;
    if ii > 20
        count = convergence(J,ii,count,error);
    end
    if  count > count_lim
        break;
    end
    
end
fprintf('Converged at %d\n',ii);
% figure;plot(1:ii,J(1:ii));

end





function shadow = findShadow(f,n_poly)
% Estimate the shadow by fitting legendre basis functions
B = LegendreBasis2D_vectorized(f,n_poly);
% Surface approximation by Legendre coeffs
P = B';
f_rep = repmat(f(:)',size(P,1),1);
P_f = P.*f_rep;
A = sum(P_f,2);
shadow = B*A;
shadow = reshape(shadow,size(f));

end




%LEGENDREBASIS compute K shifted legendre basis for the vectorized image
% j-th col of B contains P_j(x)
function [B] = LegendreBasis2D_vectorized(Img,k)

[Nr,Nc] = size(Img);
N = length(Img(:));     % Vectorized image

B = zeros(N,(k+1).^2);
[~,B_r] = legendre_1D(Nr,k);
[~,B_c] = legendre_1D(Nc,k);

ind = 0;
for ii = 1 : k+1
    for jj = 1 : k+1
        ind = ind+1;
        row_basis = B_r(:,ii);
        col_basis = B_c(:,jj);
        outer_prod = row_basis*col_basis';  % same size as the image
        B(:,ind) = outer_prod(:);
        
    end
end

end

function [B,orthonormal_B] = poly_1D(N,k)

X = 0:N-1;
p0 = ones(1,N);

B = zeros(N,k+1);
B(:,1) = 1;
orthonormal_B = B;


for ii = 2 : k+1
    B(:,ii) = B(:,ii-1).*X';
    orthonormal_B(:,ii) = B(:,ii)/norm(B(:,ii));
end

end


function [B,orthonormal_B] = legendre_1D(N,k)

X = -1:2/(N-1):1;
p0 = ones(1,N);

B = zeros(N,k+1);
orthonormal_B = B;
B(:,1) = p0';
orthonormal_B(:,1) = B(:,1)/norm(B(:,1));

for ii = 2 : k+1
    Pn = 0;
    n = ii-1;   % degree
    for k = 0 : n
        Pn = Pn +  (nchoosek(n,k)^2)*(((X-1).^(n-k)).*(X+1).^k);
    end
    B(:,ii) = Pn'/(2)^n;
    orthonormal_B(:,ii) = B(:,ii)/norm(B(:,ii));
end

end


function count = convergence(J,its,count,error)
m = J(its);
n = J(its-1);
if abs(m-n) <= error
    count = count + 1;
else
    count = 0;
end
end

