function [SetCells] = getCellParts(neuronBW,dirtSize)
%INPUTS
%neuronBW - black and white image of the neuron from function neuronBW
%dirtSize - number of pixels to consider as dirt and not part of neuron
%
%OUTPUT
%SetCells - returns images of each cell part
%  cell of {number of parts, 2}
%  stores the image of each part in the first col and the x,y location of
%  the left top corner in the second col.


[partsCell, num] = bwlabel(neuronBW);
ind = cell(num,1);
for n = 1:num
    ind{n} = find(partsCell == n); % search for the n-th conn comp.
    pts(n) = length(ind{n});
end

f = find(pts > dirtSize);

CleanCell = zeros(size(neuronBW));
% discard those conn comps 
for m = 1:length(f)
    CleanCell(ind{f(m)}) = 1;
end

% CleanCell = bwmorph(CleanCell,'erode',1);

% get cell parts
[curLabel, num] = bwlabel(CleanCell);

ind2 = cell(num,1);
for n = 1:num
    ind2{n} = find(curLabel==n);
    [rcell{n} ccell{n}] = ind2sub(size(neuronBW),ind2{n});
end

%make little images
SetCells = cell(num,2);
for m = 1:num
    minr = min(rcell{m});
    minc = min(ccell{m});
    tempI = zeros(max(rcell{m})-minr+1,max(ccell{m})-minc+1);
    for p = 1:length(rcell{m})
       tempI(rcell{m}(p)-minr+1,ccell{m}(p)-minc+1) = 1;
    end
    SetCells{m,1} = tempI;
    SetCells{m,2} = [minc, minr];
end
