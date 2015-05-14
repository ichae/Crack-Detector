function [phi0] = initRegionOtsu(I,level,area)
%INITREGIONOTSU : Initialize the level set function using otsu's method
bwI = I > level*graythresh(I);
% bwI = ~bwI; % cracks are darker
bwI = bwareaopen(bwI,area);
phi0 = TuFF.computeSDF(bwI);

end

