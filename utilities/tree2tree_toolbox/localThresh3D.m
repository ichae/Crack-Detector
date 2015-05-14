function [bwI] = localThresh3D(Img,vesselImg)
%LOCALTHRESH3D Perform a local thresholding based on the image statistics


[row col depth] = size(vesselImg);
Img = Img(:,:,2:depth+1);

p = 6;

I = Img;
% I = vesselImg;
T = I;
num_iter = 200;

obtGradient = findGradient(I);
obtGradient = obtGradient.^p;

g = (obtGradient)/(max(obtGradient(:)));

lambda = 0;
prod = zeros(size(T));
gradT = zeros(size(T));

for iter = 1 : num_iter
 
        prod = (lambda*vesselImg + (1-lambda)*g).*(I-T);

        gradT = findGradient(T);
        gradT = gradT.^2;
%         gradT = gradT/max(gradT(:)); 
        E1 = 0.5*sum(prod(:));
        E2 = 0.5*sum(gradT(:));
        alpha = E2/sqrt(E1^2+E2^2);
        Laplacian = find3DLaplacian(T);
        if mod(iter,20) == 0
            disp(iter);
        end
        T = T+sqrt(1-alpha*alpha)*prod + alpha*(Laplacian);
        alpha
end

fact = .5; % 0.15 USING I = vesselImg
bwI = I > fact*T ;


end


function [grad] = findGradient(Img)

[dX dY dZ] = meshgrid([-1 0 1]);

Gx = (1/27)*imfilter(Img,dX,'conv','same');
Gy = (1/27)*imfilter(Img,dY,'conv','same');
Gz = (1/27)*imfilter(Img,dZ,'conv','same');
grad = sqrt(Gx.^2+Gy.^2+Gz.^2);


end

function [laplacian3D] = find3DLaplacian(Img)

[dX dY dZ] = meshgrid([-1 0 1]);

Gx = (1/27)*imfilter(Img,dX,'conv','same');
Gy = (1/27)*imfilter(Img,dY,'conv','same');
Gz = (1/27)*imfilter(Img,dZ,'conv','same');

Gxx = (1/27)*imfilter(Gx,dX,'conv','same');
Gyy = (1/27)*imfilter(Gy,dY,'conv','same');
Gzz = (1/27)*imfilter(Gz,dZ,'conv','same');

laplacian3D = Gxx+Gyy+Gzz;
laplacian3D = laplacian3D/max(laplacian3D(:));
end