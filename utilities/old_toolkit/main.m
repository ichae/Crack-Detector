
%% Basic Level Set Implementation
clc
clear all
close all

currFolder = pwd;
[FileName,PathName] = uigetfile('./*.*','Input Image');
fname = [PathName FileName ];
Img = imread(fname);

Img = imresize(Img,[200 200]);
Img = padarray(Img,[20 20],'replicate');

if (length(size(Img)) == 3)
    Img = rgb2gray(Img);
end

h1 = figure(1); imshow(Img); hold all;

res = 0.5;
[init_pts] = initialize(Img,res);

plot(init_pts.y,init_pts.x,'-r');

%% initialize the 

surface_map = zeros(size(Img));
surface_map = initSurface(surface_map,init_pts);

%% 
iter = 50;
tau = 0.5;

phi = surface_map;
F = 1;  % force

for i = 1 : iter
    [gy gx] = gradient(phi);
    grad_phi = sqrt(gx.^2 + gy.^2);
    F = 1./(1+(grad_phi).^2);
%     F = 1;
    new_phi = phi + (tau)*grad_phi + tau*F;
    if mod(i,10) == 0
        contour_pts = findZeroCrossing(new_phi);
        p.x = contour_pts(:,1); p.y = contour_pts(:,2);
        new_phi = initSurface(new_phi,p);
        
        plot(contour_pts(:,2),contour_pts(:,1),'-'); drawnow;hold all;
        
        
    end
%     imagesc(new_phi);title(num2str(i));hold on;drawnow;
    disp(i)
    phi = new_phi;
    
end