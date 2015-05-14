% Code by Suvadip
% Neuron Segmentation
%% Initiate inside the structure and then evolve

clc
clear all
close all

addutilities();
pad_size            = 0;
resize              = 1;
multiple_init       = 1;
Img = readImage(pad_size,resize);
% Img = medfilt2(Img,[3 3]);
[nrow, ncol] = size(Img);
% imshow(Img);

%% Shadow removal
se = strel('disk',3);
clImg = imclose(Img,se);
sigma = 3;
sz = 7*sigma;
f = fspecial('gaussian',sz,sigma);
clImg = conv2(clImg,f,'same');
% clImg = medfilt2(clImg,[20 20]);
crackImg = clImg - Img;
se2 = strel('disk',1);
% crackImg = imclose(crackImg,se2);
crackImg = medfilt2(crackImg,[3 3]);
figure(1); 
subplot(1,2,1); imshow(clImg);
subplot(1,2,2); imshow(crackImg);




%% Surface evolution -- Parameter setting

initial_inside = 1;

smin = 3;
smax = 3;
[F1,Hessian, enhImg] = obtainVesselTangentFieldAndHessian(crackImg,smin,smax);
enhImg = (enhImg - min(enhImg(:)))/(max(enhImg(:))-min(enhImg(:)));

%%
figure(1);imshow(enhImg);
h= imrect;
mask = createMask(h);
M = enhImg.*mask;
% close(1);

se = strel('disk',3);
M_close = imclose(M,se);
M_open  = imopen(M,se);
M_oc    = imclose(M_open,se);
M = M_open;
M = (M-min(M(:)))/(max(M(:))-min(M(:)));

%% Surface initialization

if multiple_init
    binary_img = im2bw(M,0.2*graythresh(M));
    open_sz = 20;
    binary_img = bwareaopen(binary_img,open_sz,8);
    imshow(binary_img);
    surface_map =  double(bwdist(~binary_img) - bwdist(binary_img));
    
else
    prop.Img = Img;
    prop.n_regions = 3;
    prop.mask_type = 'single';
    prop.img_magnify = 100;
    init_mask = initRegions(prop);
    surface_map = computeSDF(init_mask);
end


%% Evolve Surface

figure(2);

epsilon                         = 5e-15;
phi                             = double(surface_map);
p_phi                           = phi;
num_iter                        = 800;
reinit_interval                 = 1;
band_rad                        = 2;
upper_bound                     = 30;
stored_max_change               = [];
change_deviation                = Inf;
mu                              = 0.0001;

se = strel('disk',3);
Dirac     = @(x,e) ((e/pi)./(e^2.+ x.^2));
% Dirac     = @(x,e) ((1/2/e)*(1+cos(pi*x/e)).*((x<=e) & (x>=-e)));
Heaviside = @(y,e) (0.5*(1+(2/pi)*atan(y/e)));


count = 0;
its = 0;
stop = 0;
for iter = 1 : num_iter
    
    [phi_x, phi_y]          = gradient(phi);           
    grad_phi_mag            = sqrt(phi_x.^2+phi_y.^2)+ epsilon;
    grad_phi(:,:,1)         = phi_x;
    grad_phi(:,:,2)         = phi_y;
    [nxx,dummy]             = gradient(phi_x./(grad_phi_mag));
    [dummy,nyy]             = gradient(phi_y./(grad_phi_mag));
    delta_phi               = Dirac(phi,band_rad);
%     delta_phi = grad_phi_mag;
    
    change =  (M.*delta_phi + mu*delta_phi.*(nxx+nyy)).*mask ;
                  
    t_step  = 0.7/max(change(:));
    p_phi = phi;
    phi = phi + t_step*change;
    phi = NeumanBoundCond(phi);
%     phi = imopen(phi,se);
    phi = reinitialization(phi,0.2);
    
    count = convergence(p_phi>0,phi>0,0.002,count);
    % count how many succesive times we have attained convergence, reduce local minima
    if count > 100
        stop = 1;
        fprintf('Algorithm converged, iteration=%d \n',iter);
    end
    if mod(iter,30) == 0
%         imagesc(phi);   drawnow;
        showCurveAndPhi(phi,enhImg,'r'); 
        drawnow;
        disp(iter)
    end

    if stop
        break;
    end
end


%%
binary_segment = mask.*(phi>0);
preserve_one_cc  = 0;
binary_segment = bwareaopen(binary_segment,1);
se = strel('disk',2);
binary_segment = imclose(binary_segment,se);
if preserve_one_cc
    CC = bwconncomp(binary_segment);      % preserve the largest component only
    numCC = CC.NumObjects;
    
    t = 0;
    for i = 1 : numCC
        sz = length(CC.PixelIdxList{i});
        if sz > t
            t = sz ;
        end
    end
    binary_segment = bwareaopen(binary_segment,t-2,8);
end


% Plot the skeleton

skel = bwmorph(binary_segment,'skel',Inf);
skel = bwmorph(skel,'spur',10);
skel_graph = skel2graph(skel);
skeleton_graph = sparse(skel_graph.mat);
node_coord     = skel_graph.node_coord;
node_id        = skel_graph.node_id; 
skeleton_graph = graphminspantree(skeleton_graph,1,'Method','Kruskal');
[R,C] = ind2sub(size(skel),node_coord);
imshow(Img);hold on;
% gplot(skeleton_graph,[C R],'--r');
plot(C,R,'*b','MarkerSize',3);
hold off;






