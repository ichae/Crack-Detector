
clear all
close all

end_frame=100;
skip=1;

I_dummy=rgb2gray(imread('JJC-12_0002_000_small.jpg'));
[r c]=size(I_dummy);
proj=zeros(r,c);

for I=1:skip:end_frame,
    
    filename=sprintf('JJC-12_0002_%0.3d_small.jpg',I);
    I_dummy=imread(filename);
    proj=proj+I_dummy;
    
end

maxP=max(proj(:));
minP=min(proj(:));

P=((proj-minP)./(maxP-minP))*255;

P=uint(P);

imshow(P);

