function [block , start , finish , xL,yL,zL] = createSubImg3D(I,offset,zoffset,xmin,ymin,zmin,xmax,ymax,zmax,Xact,Yact,Zact)
%CREATESUBIMG3D Create a subImage

[row col depth] = size(I);

% zoffset = offset;

xL = xmin - offset; 
if xL < 1
    xL = 1;
end

yL = ymin - offset ;
if yL < 1
    yL = 1;
end

zL = zmin - zoffset ;
if zL < 1
    zL = 1;
end

xR = xmax + offset; 
if xR > row
    xR = row;
end

yR = ymax + offset ;
if yR > col
    yR = col;
end

zR = zmax + zoffset ;
if zR > depth
    zR = depth;
end

% xmin = min(xL,xR);xmax = max(xL,xR);
% ymin
% ymax

block = I(xL:xR,yL:yR,zL:zR);

[row col depth] = size(block);

start(1) = Xact(1) - xL;
if start(1) < 1
    start(1) = 1;
elseif start(1)>row
    start(1) = row;
end

start(2) = Yact(1) - yL;
if start(2) < 1
    start(2) = 1;
elseif start(2)>col
    start(2) = col;
end

start(3) = Zact(1) - zL;
if start(3) < 1
    start(3) = 1;
elseif start(3)>depth
    start(3)=depth;
end

finish(1) = Xact(2) - xL;
if finish(1)> row
    finish(1) = row;
elseif finish(1)<1
    finish(1) = 1;
end

finish(2) = Yact(2) - yL;
if finish(2)> col
    finish(2) = col;
elseif finish(2)<1
    finish(2)=1;
end

finish(3) = Zact(2) - zL;
if finish(3)> depth
    finish(3) = depth;
elseif finish(3)<1
    finish(3) = 1;
end



end

