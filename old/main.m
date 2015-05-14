% Create the template
clc
clear all
close all
inImg = imread('crack2.jpg');
dim = ndims(inImg);

if(dim == 3)
    %Input is a color image
%     inImg = imnoise(inImg,'gaussian',0,0.01);
    inImg = rgb2gray(inImg);
end
inImg = im2double(imresize(inImg,1.5));
I = medfilt2(inImg);
I = inImg;
[nrow,ncol] = size(I);
figure(1);imshow(I);


%%
%------ Tunables -------------
drawtemplate        = 0;
sgn                 = 1;  % -1 if the vessel is bright
sigma               = 3;
d                   = 1.2*sigma;
% d                   = -1 ;        % signifies use only local detection
all_theta           = - 90:10:90;
psi                 = -20:5:20;
%------ Fixed -----------------
l_theta             = length(all_theta);
l_psi               = length(psi);
ind                 = 0;
siz     = round(8*sigma);
[X,Y]   = meshgrid(-siz:siz);
G       = exp(-(X.^2+Y.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
Gxx     = (X.^2-sigma^2).*G/(sigma^4);
Gyy     = (Y.^2-sigma^2).*G/(sigma^4);
Gxy     = (X.*Y).*G/(sigma^4);


for ii = 1 : l_theta
    rot = all_theta(ii);
    substack = [];
    for jj = 1:l_psi
        ind = ind+1;
        orientation = psi(jj)+rot;
        if d == -1
            X_d = 0;
            Y_d = 0;
            G_d = 0;
            Gxx_d = 0;
            Gyy_d = 0;
            Gxy_d = 0;
        else
            X_d = X - d*sind(orientation);
            Y_d = Y + d*cosd(orientation);
            G_d = exp(-(X_d.^2+Y_d.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
            Gxx_d = (X_d.^2-sigma^2).*G_d/(sigma^4);
            Gyy_d = (Y_d.^2-sigma^2).*G_d/(sigma^4);
            Gxy_d = (X_d.*Y_d).*G_d/(sigma^4);
        end
        
        R_theta = Gxx*(cosd(rot))^2 + Gyy*(sind(rot))^2 + Gxy*(sind(2*rot))+...
            +(Gxx_d)*(cosd(orientation))^2+(Gyy_d)*(sind(orientation))^2 + (Gxy_d)*(sind(2*orientation));
        I_theta = imfilter(I,sgn*sigma^2*R_theta,'same','replicate');
        substack(:,:,jj) = I_theta;
        if drawtemplate
            figure(2);
            subtightplot(l_theta,l_psi,ind,0.01,0.01,0.01);imagesc(R_theta);
            drawnow; colormap('gray'); axis off;
        end
    end
    stack(:,:,ii) = max(I_theta,[],3);
end

%%
% Find the total response
thresh_red          = 1;
resp = max(stack,[],3);
% resp = (resp-min(resp(:)))/(max(resp(:))-min(resp(:)));
figure(3);imshow(resp,'InitialMagnification',100);
bwI = resp > thresh_red*graythresh(resp);
bwI = bwareaopen(bwI,400);
se = strel('disk',1);
bwI = imclose(bwI,se);
figure(5);imshow(bwI,'InitialMagnification',100);

%% Plot the skeleton
skel = bwmorph(bwI,'skel',Inf);
skel = bwmorph(skel,'spur',5);
[R,C] = find(skel);
figure(6); imshow(I,[],'InitialMagnification',100); hold on;
plot(C,R,'.g','MarkerSize',5); 
hold off;



