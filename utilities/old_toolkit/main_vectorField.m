%% Initiate inside the structure and then evolve

clc
clear all
close all

pad_size = 0;
resize = 1;
Img = readImage(pad_size,resize);

% Img = Img > 0.5;
h1 = figure(1); imshow(Img); hold all;
res = 0.5;
[init_pts] = initialize(Img,res);
plot(init_pts.y,init_pts.x,'-r');

%% Initialize the surface
surface_map = zeros(size(Img));
surface_map = -initSurface(surface_map,init_pts);

[nrow ncol] = size(Img);

%% Surface evolution -- Parameter setting

sigma = 2;
h = fspecial('gaussian',ceil(12*sigma),sigma);
smooth_Img = conv2(Img,h,'same');

[gx gy] = gradient(smooth_Img);
grad_Img = sqrt(gx.^2+gy.^2);
% rho = 5;
% outward_force = exp(-grad_Img/rho);

flux_divergence = (divergence(gx,gy));
t = find(flux_divergence>0);
flux_divergence(t) = 0;
flux_divergence = abs(flux_divergence);
outward_force = flux_divergence;
rho = 3;
% outward_force = exp(-flux_divergence/rho);

% Surface evolution -- Optimization

h2 = figure(2); imshow(Img); hold all;
phi = surface_map;
new_phi = phi;
t_step = 2;
num_iter = 500;

mu = 0;
initial_inside = 1;

F = outward_force;
% F = outward_force.*Img;
F = (F - min(F(:)))/(max(F(:))-min(F(:)));
p.x = [];
p.y = [];
for iter = 1 : num_iter

    [phi_x phi_y] = gradient(phi);
    grad_phi = sqrt(phi_x.^2+phi_y.^2);
    change = t_step*grad_phi.*F ;
    phi = phi + change;
    if mod(iter,5) == 0
        [contour_pts,phi] = findZeroCrossing(phi,p,1);
        p.x = contour_pts(:,1); p.y = contour_pts(:,2);
%         plot(p.y,p.x,'-'); drawnow;hold all; 
        imshow(phi > 0); hold on; drawnow;
        disp(iter);
    end
    
end



