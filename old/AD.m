
clc
clear all
close all
inImg = imread('sidewalk_iphone.jpg');
dim = ndims(inImg);
mag = 100;
rgbImg = inImg;
if(dim == 3)
    colImg = inImg;
    inImg = rgb2gray(inImg);
end
inImg = im2double(inImg);
inImg = imresize(inImg,0.4);

% [nrow, ncol] = size(I);
% figure(1);imshow(I,'InitialMagnification',mag);

%% Perform Anisotropic diffusion

iter = 50;
k = 0.02;
u = inImg;
% u = imnoise(u,'gaussian',0.1);
figure(1);imshow(u);
for jj =  1 : iter
    [gx,gy] = gradient(u);
%     [gxx,~] = gradient(gx);
%     [~,gyy] = gradient(gy);
    g = sqrt(gx.^2+gy.^2);
%     g = g./(max(g(:)));
    c = 1./(1+(g/k).^2);
    c = c./(max(c(:)));
%     c = 1;
    [gxx,~] = gradient(c.*gx);
    [~,gyy] = gradient(c.*gy);
    term = gxx+gyy;
%     term = term/(max(abs(term(:)))+0.1);
    dt = 0.1/(max(abs(term(:)))+0.1);
    u = u + dt*term;
    u = NeumanBoundCond(u);
    imshow(u);drawnow;
    
    
end
%%
s = 64;
sz = 6*s;
f = fspecial('gaussian',sz,s);
I_f = imfilter(inImg,f,'same','replicate');
I = I_f-u;
subplot(1,3,1); imshow(inImg);
subplot(1,3,2); imshow(I_f);
subplot(1,3,3); imshow(I);


