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

sigma = 4;
h = fspecial('gaussian',ceil(12*sigma),sigma);
smooth_Img = conv2(Img,h,'same');

[gx gy] = gradient(smooth_Img);
grad_Img = sqrt(gx.^2+gy.^2);
rho = 5;
outward_force = exp(-grad_Img/rho);

%% Surface evolution -- Surface evolution

h2 = figure(2); imshow(Img); hold all;
phi = surface_map;
new_phi = phi;
t_step =3;
num_iter = 500;

mu = 0;
initial_inside = 1;

F = outward_force.*Img;
F = (F - min(F(:)))/(max(F(:))-min(F(:)));
for iter = 1 : num_iter

    laplacian_phi = del2(phi);
    [phi_x phi_y] = gradient(phi);
    grad_phi = sqrt(phi_x.^2+phi_y.^2);
    
    divergence_term = divergence(phi_x./(grad_phi+0.1),phi_y./(grad_phi+0.1));
    
    change = t_step*grad_phi.*F + t_step*mu*(laplacian_phi + divergence_term);
    
    new_phi = phi + change;
    
    if mod(iter,1) == 0
        [contour_pts,new_phi] = findZeroCrossing(new_phi,1);
        p.x = contour_pts(:,1); p.y = contour_pts(:,2);
        plot(p.y,p.x,'-'); drawnow;hold all; 
        disp(iter);
        
    end
    
    phi = new_phi;
end



