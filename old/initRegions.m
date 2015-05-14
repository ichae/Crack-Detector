function [bwI] = initRegions(prop)
%CREATEMASK Initialize the mask for initialization
I               = prop.Img;
N               = prop.n_regions;
Type            = prop.mask_type;
mag             = prop.img_magnify;

mask = false(size(I));
figure(1);imshow(I,'InitialMagnification',mag); hold on;
if strcmp(Type,'rect')
    for ii = 1 : N
        m(ii) = imrect;
        mask = mask + createMask(m(ii));
    end
elseif  strcmp(Type,'ellipse')      % default type
    for ii = 1 : N
        m(ii) = imellipse;
        mask = mask + createMask(m(ii));
    end
elseif strcmp (Type,'multiball')

    v = maskcircle2(I,'whole');
    mask = v(:,:,1);
    
elseif strcmp (Type,'single')
    
    init_rad = 5;
    mask       = zeros(size(I));
    hold on;
    for ii = 1 : N
        [init_pts] = initialize(I);
        mask = mask+ initSurface(I,init_pts , init_rad);
    end
    hold off;
else    % Default is rectangles
    for ii = 1 : N
        m(ii) = imrect;
        mask = mask + createMask(m(ii));
    end
end

close(1);
bwI = mask;

end

function m = maskcircle2(I,type)
% auto pick a circular mask for image I 
% built-in mask creation function
% Input: I   : input image     
%        type: mask shape keywords
% Output: m  : mask image

% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  

if size(I,3)~=3
    temp = double(I(:,:,1));
else
    temp = double(rgb2gray(I));
end

h = [0 1 0; 1 -4 1; 0 1 0];
T = conv2(temp,h);
T(1,:) = 0;
T(end,:) = 0;
T(:,1) = 0;
T(:,end) = 0;

thre = max(max(abs(T)))*.5;
idx = find(abs(T) > thre);
[cx,cy] = ind2sub(size(T),idx);
cx = round(mean(cx));
cy = round(mean(cy));

[x,y] = meshgrid(1:min(size(temp,1),size(temp,2)));

m = zeros(size(temp));
[p,q] = size(temp);

switch lower (type)
    case 'small'
        r = 10;
        n = zeros(size(x));
        n((x-cx).^2+(y-cy).^2<r.^2) = 1;
        m(1:size(n,1),1:size(n,2)) = n;
        %m((x-cx).^2+(y-cy).^2<r.^2) = 1;
    case 'medium'
        r = min(min(cx,p-cx),min(cy,q-cy));
        r = max(2/3*r,25);
        n = zeros(size(x));
        n((x-cx).^2+(y-cy).^2<r.^2) = 1;
        m(1:size(n,1),1:size(n,2)) = n;
        %m((x-cx).^2+(y-cy).^2<r.^2) = 1;
    case 'large'
        r = min(min(cx,p-cx),min(cy,q-cy));
        r = max(2/3*r,60);
        n = zeros(size(x));
        n((x-cx).^2+(y-cy).^2<r.^2) = 1;
        m(1:size(n,1),1:size(n,2)) = n;
        %m((x-cx).^2+(y-cy).^2<r.^2) = 1;
    case 'whole'
        r = 5;
        m = zeros(round(ceil(max(p,q)/2/(r+1))*3*(r+1)));
        siz = size(m,1);
        sx = round(siz/2);       
        i = 1:round(siz/2/(r+1));
        j = 1:round(0.9*siz/2/(r+1));
        j = j-round(median(j)); 
        m(sx+2*j*(r+1),(2*i-1)*(r+1)) = 1;
        se = strel('disk',r);
        m = imdilate(m,se);
        m = m(round(siz/2-p/2-6):round(siz/2-p/2-6)+p-1,round(siz/2-q/2-6):round(siz/2-q/2-6)+q-1); 
end
        tem(:,:,1) = m;
        M = padarray(m,[floor(2/3*r),floor(2/3*r)],0,'post');
        tem(:,:,2) = M(floor(2/3*r)+1:end,floor(2/3*r)+1:end);
        m = tem; 
end