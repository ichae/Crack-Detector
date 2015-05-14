clc
clear all
close all
folder = './';
imname = 'crack2';
extn = '.jpg';
inImg = imread(strcat(folder,imname,extn));
dim = ndims(inImg);
mag = 100;

inImg = imresize(inImg,2);
if(dim == 3)
%     inImg = imnoise(inImg,'gaussian',0,0.02);
    colImg = inImg;
    inImg = rgb2gray(inImg);
end
inImg = im2double(inImg);
[nrow, ncol] = size(inImg);

%% Remove shadow + clutter
close all;

inImg = (inImg-min(inImg(:)))/(max(inImg(:))-min(inImg(:)));
param.Img = (inImg);          % signal to denoise+deshadow
param.type = 'minimax';     % gaussian/legendre/none(only denoising)
param.n_poly = 4;          % #poly for legendre curve fitting  
param.sigma =  15;          % sigma for gaussian approximation
param.n_iter = 0;             
param.lambda = 0.6;           % smoothness coefficient
param.error = 0.005;         % convergence criteria
param.cnvg_cnt = 5;        % steady state convg. criteria
param.intrvl = 50;          % display interval

[smoothImg,shadow] = shadowRemove2(param);
%% Plot
figure(1); 
subplot(1,2,1); imshow((inImg),[]);
subplot(1,2,2); imshow(smoothImg);

