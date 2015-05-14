function [histeq_Img] = imhist3D(Img)
%IMHIST3D Compute histogram equalization of 3D images
% Img                   : Input image, type double
% histeq_Img            : image with equalized histogram

[nrow,ncol,nslice]  = size(Img)     ;
Img                 = round(Img*255);        % convert to 8-bit image
nbins               = 256           ;
bins                = 0:nbins-1     ;

[histogram]         = hist(Img(:), bins)  ;        % obtain the histogram
histogram = histogram/length(Img(:));        % normalized histogram
CDF                 = cumsum(histogram);

histeq_Img          = zeros(size(Img));

for ii = 1 : nrow
    for jj = 1 : ncol
        for kk = 1: nslice
            histeq_Img(ii,jj,kk) = CDF(Img(ii,jj,kk)+1);
        end
    end
end

end

