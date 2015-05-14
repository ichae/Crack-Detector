
function segI = binarize3D(origImg,Img,tholdLevelReduce,openSz,minObjSize,bndCut,threshOption)
% BINARIZE the 3D enhanced Image
%       USAGE: segI = binarize3D(Img,tholdLevelReduce,openSz,minObjSize,bndCut)
%       Img: Enhanced Image
%       tholdLevelReduce: Fraction of the threshold valu obtained from Otsu
%       openSz: Size of the structuring element to to open
%       minObjSize: minimum size of a connected component to keep
%       bndCut: Cut away the boundary to get rid of boundary effects


mn = min(Img(:));
mx = max(Img(:));
Img = (Img - mn)/(mx-mn);

if threshOption == 1
    level = tholdLevelReduce*graythresh(Img(:));
    tImg = zeros(size(Img));
    tImg(Img >= level) = 1;

elseif threshOption == 2
    tImg = localThresh3D_new(origImg,Img,tholdLevelReduce);
%     tImg = origImg > 1.2*T;
    
elseif threshOption == 3
    tImg = textureLocalThresh(origImg,Img,tholdLevelReduce);
end

    
    % imopen for dust
    
    disp('imopen3D!');
    segI = imopen(tImg,ones(openSz,openSz,openSz));
%      disp('imoclose3D!');
    segI = imclose(tImg,ones(3,3,3));
    % Remove very small connected components
    
    [L,NUM] = bwlabeln(segI);
    objSize = zeros(1,NUM);
    
    for I=1:NUM,
        
        dummy = any(L == I);
        objSize(I)=sum(dummy(:));
        
    end
    
    % remove the objects smaller than minimum mentioned size
    for I=1:NUM,
        
        if objSize(I) < minObjSize,
            
            segI(L==I) = 0;
            
        end
        
    end





% Get rid of boundary effects

segIdummy = zeros(size(segI));
segIdummy(bndCut+1:size(segI,1)-bndCut+1,bndCut+1:size(segI,2)-bndCut+1,bndCut+1:size(segI,3)-bndCut+1)=...
segI(bndCut+1:size(segI,1)-bndCut+1,bndCut+1:size(segI,2)-bndCut+1,bndCut+1:size(segI,3)-bndCut+1);
segI=segIdummy;