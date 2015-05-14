
function Img = read3DStack(sliceReduce,stackSkip)
% READ the 3D stacks interactively
%       USAGE: I = read3DStack(sliceReduce,stackSkip);
%       I: the 3D image after Full Scale Contrast Stretch 
%       sliceReduce: reduce the image size , stackSkip: number of intermediate
%       stacks to skip

dname = uigetdir('','Input the image folder');
D = dir(strcat(dname,'\*.tif'));
buffer = double(imread(strcat(dname,'\',D(1).name)));
buffer = imresize(buffer,sliceReduce);
[maxr maxc] = size(buffer);
Img = zeros(maxr,maxc,length(1:stackSkip:length(D)));
count=0;

for I=1:stackSkip:length(D),
    count=count+1;
    buffer = im2double(imread(strcat(dname,'\',D(I).name)));
    buffer = imresize(buffer,sliceReduce);
    Img(:,:,count) = buffer;
    
end

if sliceReduce ~=1 && stackSkip ~= 1
    Img = (Img-min(Img(:)))/(max(Img(:))-min(Img(:)));
end