function [enhI] = vesselEnhanceLDE(param)
%VESSELENHANCELDE Enhance the vessels via LDE

vessel_type = param.vessel_type;        %-- 'bright' or 'dark' vessel
I     = param.Img;                      %-- original signal
theta = param.detector_orientation;     %-- orientation space for detectors    
psi   = param.predictor_orientation;    %-- orientation space for predictors
sigma = param.scale;                    %-- scale space
k     = param.offset_factor;            %-- d = k*sigma

if strcmp(vessel_type,'dark')
    sgn = 1;
else
    sgn = -1;
end
[nr,nc] = size(I);
stack = zeros(nr,nc,length(sigma));

for ii = 1 : length(sigma)
    d = k*sigma(ii);
    stack(:,:,ii) = LDE(I,theta,psi,d,sigma(ii),sgn);
end

enhI = max(stack,[],3);

end



function [resp] = LDE(I,all_theta,psi,d,sigma,sgn)
%LDE Perform local directional evidence filtering
%   all_theta -- orientation for detection filter
%   psi       -- orientations of evidence filter, for each value of theta
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
        I_theta = (I_d + 1*I_b + 1*I_f) ;                         % superposition
% %                I_theta = 1*max(cat(3,I_d,I_f,I_b),[],3) ;     % max oriented response
%         t = find(I_theta>1);
%         I_theta(t) = 1;
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

resp = max(stack,[],3);

end

