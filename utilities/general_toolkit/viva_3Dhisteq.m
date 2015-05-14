function [finalI] = viva_3Dhisteq(Img)
%3DHISTEQ Perform the histogram equalization of a 3D Image
%       Img: A 3D grayscale image; finalI: the histogram equalized image
%       Author: Suvadip Mukherjee(sm5vp@virginia.edu)

[row col depth] = size(Img);
num_pixels = row*col*depth;

Img = round(Img.*255); % convert the gray values between 0 and 255

histogram = zeros(256,1);

vector_Img = Img(:);

% Create the histogram
for i = 1:length(vector_Img)
   histogram(vector_Img(i)+1) = histogram(vector_Img(i)+1)+1; 
end

histogram = histogram/255;

% Form The cumulative Histogram
cum_hist = zeros(256,1);

cum_hist(1) = histogram(1);

for i = 2:256
    cum_hist(i) = cum_hist(i-1)+histogram(i);
end

% create the histogram Image
finalI  = zeros(row,col,depth);

for i = 1:row
    for j = 1 : col
        for k = 1 : depth
            finalI(i,j,k) = cum_hist(Img(i,j,k)+1);
        end
    end
end

% finalI = finalI/num_pixels;
finalI = im2double(finalI);

finalI = (finalI-(min(finalI(:))))/(max(finalI(:))-min(finalI(:)));

end

