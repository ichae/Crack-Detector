function [Img] = readImage(pad_sz,do_resize)
%READIMAGE Read a 2D image
% currFolder = pwd;
[FileName,PathName] = uigetfile('./*.*','Input Image');
fname = [PathName FileName ];
Img = imread(fname);

if do_resize ~= 0
    Img = imresize(Img,do_resize);
end
Img = padarray(Img,[pad_sz pad_sz],'replicate');

if (length(size(Img)) == 3)
    Img = rgb2gray(Img);
end

Img = double(Img);

Img = (Img - min(Img(:)))/(max(Img(:))-min(Img(:)));
end

