function [s] = estimateShadowMinimax(obj,deshadow_param)
%ESTIMATESHADOWMINIMAX Shadow estimation for thin structures using
%variational minimax
close all;

p = deshadow_param.power;                   %-- exponent of gradient magnitude
sigma = deshadow_param.sigma;               %-- scale for initial surface
n_iter = deshadow_param.n_iter_shadow;      %-- number of iterations
error = deshadow_param.error;               %-- convergence error criteria
count_lim = deshadow_param.cnvg_cnt;        %-- convergence count

f = obj.inImg;
sz = 6*sigma;
h1 = fspecial('gaussian',sz,sigma);
s = imfilter(f,h1,'same','replicate');

[gx,gy] = gradient(medfilt2(f));
g = (sqrt(gx.^2+gy.^2)).^p; 
g = g/max(g(:));
figure('name','Shadow Function');imshow(s);
J = zeros(n_iter,1);
count = 0;
for ii = 1 : n_iter
    s = DeShadow.NeumannBoundCond(s);
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
    imshow(s,[]);drawnow;
    
    % Check convergence
    J(ii) = sqrt(1-alpha^2)*E1 + alpha*E2;
    if ii > 20
        count = DeShadow.convergence(J,ii,count,error);
    end
    if  count > count_lim
        break;
    end
    
end
fprintf('Converged at %d\n',ii);




end

