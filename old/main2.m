% Create the template
clc
clear all
close all
inImg = imread('sidewalk_iphone.jpg');
dim = ndims(inImg);
mag = 100;
% rgbImg = inImg;
if(dim == 3)
    %Input is a color image
%     inImg = imnoise(inImg,'gaussian',0,0.02);
    colImg = inImg;
    inImg = rgb2gray(inImg);
end
inImg = im2double(inImg);
inImg = imresize(inImg,0.5);

% I = inImg;
[nrow, ncol] = size(inImg);
figure(1);imshow(inImg,'InitialMagnification',mag);

%%
s = 20;
sz = 6*s;
f = fspecial('gaussian',sz,s);
I_f = imfilter(inImg,f,'same','replicate');
% I = medfilt2(inImg,[15 15]);
I = I_f-inImg;
subplot(1,3,1); imshow(inImg);
subplot(1,3,2); imshow(I_f);
subplot(1,3,3); imshow(I);

%% Local Directional Evidence
%------ Tunables -------------

sgn                 = 1;  % -1 if the vessel is bright
sigma               = 4;
d                   = 0.8*sigma;
% d                   = -1 ;        % signifies use only local detection
all_theta           = - 90:15:90;
psi                 = -30:10:30;

resp = LDE(I,all_theta,psi,d,sigma,sgn);
% resp_MF = LDE(I,all_theta,psi,-1,sigma,sgn); % d=-1 signifies only  detector
% resp = (resp-min(resp(:)))/(max(resp(:))-min(resp(:)));
%%
figure; subplot(1,2,1); imshow(resp);
subplot(1,2,2);imshow(resp>1*graythresh(resp));

%%
[gx,gy] = gradient(resp);
g = sqrt(gx.^2+gy.^2);
imshow(g,[])
%% Frangi
%% Frangi filter
size = 6*sigma;
h = fspecial('gaussian',size, sigma);
gamma = 0.25;
if sgn < 0  % bright vessel
    I_smooth = conv2(I,h,'same');
    frangi_resp = sigma^gamma*vesselEnhanceBright(I_smooth,sigma);
else
    I_smooth = conv2(I,h,'same');
    frangi_resp = sigma^gamma*vesselEnhanceDark(I_smooth,sigma);
end

%% ---------------------- Visual Comparison of results
close all;
thresh_red          = 0.5;
area                = 50;
spur_amt            = 5;
% compareVisually(I,resp,frangi_resp,thresh_red,spur_amt,area,mag);
compareVisuallySeparate(inImg,resp,frangi_resp,thresh_red,spur_amt,area,mag);













