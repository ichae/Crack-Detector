
%% Reading the image
clc
clear all
close all

[colImg,inImg,imname, pathname] = StdIP.readImg2D(0.4);
% inImg = inImg/max(inImg(:));
% inImg = imadjust(inImg);
imshow(inImg);

%% Shadow + Clutter removal

deshadow_obj = DeShadow(inImg);
deshadow_param.n_iter_shadow    = 500;
deshadow_param.n_iter_reconst   = 500;
deshadow_param.error            = 0.01;
deshadow_param.cnvg_cnt         = 15;
deshadow_param.interval         = 50;
deshadow_param.smooth_trm       = 0.8;
%-- parameters to tune
deshadow_param.sigma            = 5;
deshadow_param.power            = 2.4;

sh3 = deshadow_obj.estimateShadowMinimax(deshadow_param);   
deshadow_param.shadow           = sh3;
smoothImg = deshadow_obj.reconstructImage(deshadow_param);
smoothImg = smoothImg./max(smoothImg(:));
close all;
%% Segmentation with TuFF: Vessel Enhancement 

% enhImg = TuFF.darkVesselEnhanceFrangi(smoothImg,2:4,2);     %-- static function of TuFF
param_LDE.vessel_type           = 'dark';
param_LDE.Img                   = smoothImg;
param_LDE.detector_orientation  = -90:15:90;
param_LDE.predictor_orientation = -30:10:30;
param_LDE.scale                 = [2];
param_LDE.offset_factor         = 0.6;      %-- d = offset_factor*sigma

enhImg = TuFF.vesselEnhanceLDE(param_LDE);                    %-- enhance vessels with LDE  
imshow(enhImg);
%% Segmentation with TuFF: Initialization
% phi0   = TuFF.initRegionManual(inImg, 'ellipse',200);       %-- Static function of TuFF
phi0 = TuFF.initRegionOtsu(enhImg,0.6,200);
imshow(phi0>0);
%% Segmentation with TuFF: Curve Evolution

param_Tuff.phi0         = phi0;
param_Tuff.enhI         = enhImg;
param_Tuff.max_iter     = 400;
param_Tuff.error        = 0.5;
param_Tuff.cnvg_cnt     = 10;
param_Tuff.interval     = 25;
param_Tuff.disp_Img     = colImg;
param_Tuff.magnify      = 200;
param_Tuff.col          = 'y';
%-- parameters to tune
param_Tuff.dt           = 0.8;
param_Tuff.pull_back    = 0.1;  % 1/0:do(not)use deflation; essential for nice shape
param_Tuff.p            = 2;    % exponent in 1/1+()^p
param_Tuff.smooth_trm   = 1.2;
param_Tuff.edge_trm     = 0.7;
param_Tuff.attr_trm     = 0.0;
param_Tuff.edge_range   = 7;
param_Tuff.attr_range   = 7;

TuFF_obj = TuFF(smoothImg);  %-- object for TuFF class
% phi = TuFF_obj.runFastEATuFF2D(param_Tuff);
phi = TuFF_obj.runEATuFF2D(param_Tuff);

%% Diplay segmentation results
im1 = colImg;
im2 = phi>0;
se_cl = strel('disk',2);
se_op = strel('disk',1);
im2 = imclose(im2,se_cl);
% im2 = imopen(im2,se_op);

% im2 = bwareaopen(im2,200);
StdIP.imageOverlay(im1,im2,200,[1 1 0]); %--- yellow = [1,1,0]


%% Save the other images



%% Delete the objects

delete(deshadow_obj);
delete(TuFF_obj);
