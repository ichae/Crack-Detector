function [SDF] = computeSDF(bwI)
%COMPUTESDF Create the signed distance function from the binary image
% inside >= 0, outside < 0

phi = bwdist(bwI)-bwdist(1-bwI)+im2double(bwI)-.5;
SDF = -phi ; 

end

