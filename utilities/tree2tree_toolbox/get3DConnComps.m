
function connComps = get3DConnComps(binImg3D)

% get conn comps

[L,NUM] = bwlabeln(binImg3D);
compPad=5;
connComps=cell(1,NUM);

for I=1:NUM,
    
    linIndex=find(L==I);
    [y,x,z]=ind2sub(size(binImg3D),linIndex);
    maxx=max(x);
    minx=min(x);
    maxy=max(y);
    miny=min(y);
    maxz=max(z);
    minz=min(z);
    zeroPointOffset=[minx-compPad,miny-compPad,minz-compPad];
    connComps{I}.Offset=zeroPointOffset;
    connComps{I}.binImg=zeros((maxy-miny+2*compPad),(maxx-minx+2*compPad),(maxz-minz+2*compPad));
    y_bar=y-zeroPointOffset(2);
    x_bar=x-zeroPointOffset(1);
    z_bar=z-zeroPointOffset(3);
    connComps{I}.binImg(sub2ind(size(connComps{I}.binImg),y_bar,x_bar,z_bar))=1;
        
end
    
