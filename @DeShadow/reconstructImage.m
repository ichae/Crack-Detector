function u = reconstructImage(obj,deshadow_param)
%RECONSTRUCTIMAGE Denoise + deshadow

shadow = deshadow_param.shadow;             %-- shadow function estimated
lambda = deshadow_param.smooth_trm;         %-- smoothness coeff
n_iter = deshadow_param.n_iter_reconst;     %-- max number of iterations
error = deshadow_param.error;               %-- convergence error criteria
count_lim = deshadow_param.cnvg_cnt;        %-- convergence count
interval = deshadow_param.interval;         %-- display interval

f = obj.inImg;
J = zeros(n_iter,1);
data_coeff = 1;
u = f;
figure('name','Shadow Removal+Denoising');
subplot(1,3,1);imshow(f);
subplot(1,3,2);imshow(u);
subplot(1,3,3);imshow(shadow);
count = 0;

for its = 1 : n_iter
    u = DeShadow.NeumannBoundCond(u);
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
        count = DeShadow.convergence(J,its,count,error);
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




