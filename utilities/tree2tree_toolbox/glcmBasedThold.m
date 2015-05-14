function [bwI] = glcmBasedThold(Img)
%GLCMBASEDTHOLD 

[row col depth] = size(Img);

numLevel = 64 ;

Img = round((numLevel-1)*Img);    % convert to 0-63 level

% hist(Img);
GLCM1 = horizontalGLCM(Img,numLevel);
GLCM2 = verticalGLCM(Img,numLevel);
GLCM3 = depthGLCM(Img,numLevel);

GLCM = (GLCM1+GLCM2+GLCM3);


% Chanda Chaudhury method

maxVals = zeros(numLevel,1);


for t = 1 : numLevel
    
   A = GLCM(1:t,1:t);
   statA = graycoprops(A);
   contrA = statA.Contrast/sum(A(:));
   B = GLCM(t+1:numLevel,t+1:numLevel);
   statB = graycoprops(B);
   contrB = statB.Contrast/sum(B(:));
   
   maxVals(t) = contrA+contrB;
end

[m,n] = max(maxVals);

param = 0.02;
thresh = param*(n-1) 
% thresh = 1;
bwI = Img > thresh;
% 
% stats = graycoprops(GLCM);
% 
% contr = stats.Contrast ;
% corr = stats.Correlation ;
% energy = stats.Energy;
% homogeneity = stats.Homogeneity;

end




function [hGLCM] = horizontalGLCM(I,numLevel)

hGLCM = zeros(numLevel,numLevel);
[row col depth] = size(I);

for i = 1 : row 
    for j = 1 : col - 1
        for k = 1 : depth
            v1 = I(i,j,k);
            v2 = I(i,j+1,k);
            hGLCM(v1+1,v2+1) = hGLCM(v1+1,v2+1)+1;
            hGLCM(v2+1,v1+1) = hGLCM(v2+1,v1+1)+1;
        end
    end
end

end

function [vGLCM] = verticalGLCM(I,numLevel)

vGLCM = zeros(numLevel,numLevel);
[row col depth] = size(I);
for i = 1 : row -1
    for j = 1 : col
        for k = 1 : depth
            v1 = I(i,j,k);
            v2 = I(i+1,j,k);
            vGLCM(v1+1,v2+1) = vGLCM(v1+1,v2+1)+1;
            vGLCM(v2+1,v1+1) = vGLCM(v2+1,v1+1)+1;
        end
    end
end

end


function [dGLCM] = depthGLCM(I,numLevel)

dGLCM = zeros(numLevel,numLevel);
[row col depth] = size(I);
for i = 1 : row 
    for j = 1 : col
        for k = 1 : depth-1
            v1 = I(i,j,k);
            v2 = I(i,j,k+1);
            dGLCM(v1+1,v2+1) = dGLCM(v1+1,v2+1)+1;
            dGLCM(v2+1,v1+1) = dGLCM(v2+1,v1+1)+1;
        end
    end
end

end










