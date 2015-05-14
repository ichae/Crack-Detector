function s = estimateShadowLegendre(obj,n_poly)
%ESTIMATESHADOWLEGENDRE Shadow estimation for thin structures using
%Legendre polynomials for modelling
%   Detailed explanation goes here
f = obj.inImg;
B = LegendreBasis2D_vectorized(f,n_poly);
% Surface approximation by Legendre coeffs
P = B';
f_rep = repmat(f(:)',size(P,1),1);
P_f = P.*f_rep;
A = sum(P_f,2);
s = B*A;
s = reshape(s,size(f));

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
