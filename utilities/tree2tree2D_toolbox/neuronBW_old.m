function [BWcell] = neuronBW(I,perCell,winSz)


I = (I-min(I(:)))/(max(I(:))-min(I(:)));
npxI = size(I,1)*size(I,2);
pixvals = linspace(0,1,256);
n = histc(I(:),0:1/256:1);
c_elements = cumsum(n)/npxI;
tau = pixvals(find(c_elements>=perCell,1,'first'));

%% get local t-hold
newI = padarray(I,floor(winSz/2),'replicate','both');
segI = zeros(size(I));

for j = 1:size(I,2);
    for i = 1:size(I,1)
       cWin = newI(i:i+winSz(1)-1,j:j+winSz(2)-1);
        level = graythresh(cWin);
        if level <= tau
            level = 1;
        else
           if I(i,j)>= level
               segI(i,j) = 1;
           else
               segI(i,j) = 0;
           end
        end
        
    end 
end
BWcell = segI;
return