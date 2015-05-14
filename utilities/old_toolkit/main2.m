% Level set with regularized distance field  
% Code based on Chunming Li et.al's paper in CVPR 2005

clc
clear all
close all

pad_size = 20;
Img = readImage(pad_size);

% Img = Img > 0.5;
h1 = figure(1); imshow(Img); hold all;
res = 0.5;
[init_pts] = initialize(Img,res);
plot(init_pts.y,init_pts.x,'-r');

%% Initialize the surface
surface_map = zeros(size(Img));
surface_map = initSurface(surface_map,init_pts);

[nrow ncol] = size(Img);
% surface_map(10:nrow-10,10:ncol-10) = -2;

%% Evolve The Surface

h2 = figure(2); imshow(Img); hold all;
num_iter = 1250;
step = .2;

phi = surface_map;

lambda_1 = 0.033;        % contribution of the surface regulizer
lambda_2 = 5;        % contribution of the zero levelset contour gradient
lambda_3 = .5;        % contribution of the gradient-volume

small_term = 1e-8;
new_phi = phi;

sigma = .5;
g = getImageTerm(Img,sigma);

for iter = 1 : num_iter
    
    laplacian_phi = del2(phi);
    [phi_x phi_y] = gradient(phi);
    grad_phi = sqrt(phi_x.^2+phi_y.^2);
    
    divergence_term = divergence(phi_x./(grad_phi+0.1),phi_y./(grad_phi+0.1));
    
    term1 = lambda_1*(laplacian_phi + divergence_term);
        
    term2 = lambda_2*(g.*(phi <= small_term).*divergence_term);
    
    term3 = lambda_3*(g.*(phi <= small_term));
    
    new_phi = phi + step*(term1+term2+term3);
%     imshow(new_phi < 0);
    
    if mod(iter,10) == 0
        
        [contour_pts,new_phi] = findZeroCrossing(new_phi);
        p.x = contour_pts(:,1); p.y = contour_pts(:,2);
        plot(contour_pts(:,2),contour_pts(:,1),'-'); drawnow;hold all;
%         new_phi = initSurface(new_phi,p);
%         imshow(new_phi < 0); drawnow; hold on;
         disp(iter);
%     mesh(new_phi);drawnow; hold on;
    end
   
    phi = new_phi;
    
end