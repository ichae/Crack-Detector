function [u,shadow] = shadowRemove_gupi( f,n_poly, n_iter,lambda , l2, error, count_lim , interval)
% Shadow removal using legendre functions
% $P_f(x)= P(x)f(x), P_{f,u}(x)= P(x)f(x)u(x) $

J = zeros(n_iter,1);
data_coeff = 1;

vect_Bases = LegendreBasis2D_vectorized(f,n_poly);
% vect_Bases = vect_Bases(:,2:end);
P = vect_Bases' ; % each row contains a basis
vect_f = f(:)';   % vectorize f into a row
rep_vect_f = repmat(vect_f,size(P,1),1);

P_f = P.*rep_vect_f;

% f = (f-min(f(:)))/(max(f(:))-min(f(:)));
% u = ones(size(f));
% s = 40;
% sz = 6*s;
% h = fspecial('gaussian',sz,s);
% u = imfilter(f,h,'same','replicate');
u = f;

figure(99);
subplot(1,3,1);imshow(f);
subplot(1,3,2);imshow(u);
subplot(1,3,3);imshow(u);
count = 0;

for its = 1 : n_iter
    u = NeumanBoundCond(u);
    [ux,uy] = gradient(u);
    g = sqrt(ux.^2+uy.^2+eps);
%     g = 1;
    [uxx,~] = gradient(ux./g);
    [~,uyy] = gradient(uy./g);
    
    
    vect_u = u(:)';
    rep_vect_u = repmat(vect_u,size(P,1),1);

    P_u = P.*rep_vect_u   ;    % P(x)u(x) 
    P_f_u = P_f.*rep_vect_u;   % P(x)f(x)u(x)
    
    L = sum(P_f_u,2);           % RHS of eqn
    K = sum(P_u,2)*sum(P_u,2)';
%     K = P_u*P_u'+ l2*eye(length(L)); % Gram matrix
%     A = K\L;                    % legendre coefficients
    A = L\K;
    A = A';
%     A_rep = repmat(A,1,size(P,2));
    shadow = P'*A;
    shadow = reshape(shadow,size(u));
%     shadow = 1;
% 
%     shadow = shadow/(max(abs(shadow(:))));
    
%     shadow = 1*(shadow-min(shadow(:)))/(max(shadow(:))-min(shadow(:)));
%         shadow = ones(size(u));


    term1 = shadow.*(f-u.*shadow);
    term2 = uxx+uyy;
    
%     term1 = term1./(max(abs(term1(:))));
%     term2 = term2./(max(abs(term2(:))));
    
    chng = data_coeff*term1+ lambda*term2;
    dt = 0.005/max(abs(chng(:)));
    u = u + dt*chng;
    
    if mod(its,interval) == 0
%         subplot(1,3,1);plot(A,'-*'); drawnow; hold on;
        subplot(1,3,2);imshow(u); drawnow;
        subplot(1,3,3);imshow(shadow); drawnow;
        disp(its);
    end
    
    % Show the error term
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
newJ = J(1:its);
fprintf('Converged at %d\n',its);
figure(100);
plot(1:its,newJ,'-');  hold on;drawnow; title('convergence');

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

