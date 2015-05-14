function [resp] = LDE(I,all_theta,psi,d,sigma,sgn)
%LDE Perform local directional evidence filtering
%   all_theta -- orientation for detection filter
%   psi       -- orientations of evidence filter, for each valur of theta
%   d         -- d = k*sigma
drawtemplate        = 0;
l_theta             = length(all_theta);
l_psi               = length(psi);
ind                 = 0;
siz     = round(6*sigma);
[X,Y]   = meshgrid(-siz:siz);
G       = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
Gxx     = (X.^2-sigma^2).*G/(sigma^4);
Gyy     = (Y.^2-sigma^2).*G/(sigma^4);
Gxy     = (X.*Y).*G/(sigma^4);
substack= [];
c1 = 0.75;
c2 = 0.25;
c3 = 0.25;

for ii = 1 : l_theta
    rot = all_theta(ii);
    substack = [];
    for jj = 1:l_psi
        ind = ind+1;
        orientation = psi(jj)+rot;
        if d == -1
            X_d = 0;
            Y_d = 0;
            G_d_f = 0;
            Gxx_d_f = 0;
            Gyy_d_f = 0;
            Gxy_d_f = 0;
            Gxx_d_b = 0;
            Gyy_d_b = 0;
            Gxy_d_b = 0;
        else
            X_d_f = X - d*sind(orientation);
            Y_d_f = Y + d*cosd(orientation);
            G_d_f = exp(-(X_d_f.^2+Y_d_f.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
            Gxx_d_f = (X_d_f.^2-sigma^2).*G_d_f/(sigma^4);
            Gyy_d_f = (Y_d_f.^2-sigma^2).*G_d_f/(sigma^4);
            Gxy_d_f = (X_d_f.*Y_d_f).*G_d_f/(sigma^4);
            
            X_d_b = X - d*sind(orientation);
            Y_d_b = Y - d*cosd(orientation);
            G_d_b = exp(-(X_d_b.^2+Y_d_b.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
            Gxx_d_b = (X_d_b.^2-sigma^2).*G_d_b/(sigma^4);
            Gyy_d_b = (Y_d_b.^2-sigma^2).*G_d_b/(sigma^4);
            Gxy_d_b = (X_d_b.*Y_d_b).*G_d_b/(sigma^4);
        end
        
        R_d = Gxx*(cosd(rot))^2 + Gyy*(sind(rot))^2 + Gxy*(sind(2*rot));
        R_b = Gxx_d_b*(cosd(orientation))^2+(Gyy_d_b)*(sind(orientation))^2 + (Gxy_d_b)*(sind(2*orientation));
        R_f = Gxx_d_f*(cosd(orientation))^2+(Gyy_d_f)*(sind(orientation))^2 + (Gxy_d_f)*(sind(2*orientation));
        
        I_d = imfilter(I,sgn*sigma^1.5*R_d,'same','replicate');   % local detector
        I_b = imfilter(I,sgn*sigma^1.5*R_b,'same','replicate');   % boosting -- backward
        I_f = imfilter(I,sgn*sigma^1.5*R_f,'same','replicate');   % boosting -- forward
        I_theta = (c1*I_d + c2*I_b + c3*I_f) ;                     % superposition

        t = find(I_theta>1);
        I_theta(t) = 1;
        t = find(I_theta<0);
        I_theta(t) = 0;
        substack(:,:,jj) = I_theta;
        if drawtemplate
            figure(2);
            subtightplot(l_theta,l_psi,ind,0.01,0.01,0.01);imagesc(R_d+R_b+R_f);
            drawnow; colormap('hot'); axis off;
        end
    end
    stack(:,:,ii) = max(substack,[],3);
end

r_mx = max(stack,[],3);
r_mn = min(stack,[],3);
% resp = (abs(r_mx)-abs(r_mn))/(0.1+abs(r_mx)+abs(r_mn));
resp = r_mx;
end

