function [bwI] = localThresh3D_new(Img,vesselImg,tholdReduce)
%LOCALTHRESH3D Perform a local thresholding based on the image statistics


[row col depth] = size(vesselImg);
Img = Img(:,:,2:depth+1);

p = 8;
old_alpha = 10;
epsilon = 2.5e-5;
fact = 1; 
lambda = 1;
tau = .4;
num_iter = 400;
startLevel = tholdReduce;

% I = Img;
I = vesselImg;
level = graythresh(I);
T = zeros(size(I));
T(:,:,:) = startLevel*level;

vesselImg = round(vesselImg*255) ;
vesselImg = double(vesselImg/255) ;

obtGradient = findGradient(I);
obtGradient = obtGradient.^p;

g = (obtGradient)/(max(obtGradient(:)));


prod = zeros(size(T));
gradT = zeros(size(T));
vesselness = vesselImg/(max(vesselImg(:)));


% figure(1); hold on;
% [X,Y] = meshgrid(1:100);
% surf1 = I(:,:,20);


% figure; hold on;
% % imagesc(I-T);
% filename = 'surfaceAnimation.gif';

for iter = 1 : num_iter
    
    prod = (lambda*vesselness + (1-lambda)*g).*(I-T);
    
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

    T = T+tau*(sqrt(1-alpha*alpha)*prod + alpha*(Laplacian));
  

% -------------------------------------------------------------------------
%     surf2 = T(:,:,20);
%     mesh(X,Y,surf1(1:100,1:100)-fact*surf2(1:100,1:100),'FaceColor','red'); drawnow;
%     
%     imagesc(I(:,:,12)-T(:,:,12)); colormap('jet');
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if iter == 1;
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
% -------------------------------------------------------------------------
%%
    error = old_alpha-alpha
    if abs(error) <= epsilon && iter > num_iter
        fprintf('\n Converged \n');
        break;
    else
        old_alpha = alpha;
    end
%     pause;
end


bwI = I > fact*T ;


end


function [grad] = findGradient(Img)

[y x z] = size(Img);
k = 1/27;

zShifted = zeros(size(Img));
zShifted(:,:,1) = Img(:,:,1);
zShifted(:,:,2:z) = Img(:,:,1:z-1);
Gz = k*(zShifted-Img);


yShifted = zeros(size(Img));
yShifted(1,:,:) = Img(1,:,:);
yShifted(2:y,:,:) = Img(1:y-1,:,:);
Gy = k*(yShifted-Img);

xShifted = zeros(size(Img));
xShifted(:,1,:) = Img(:,1,:);
xShifted(:,2:x,:) = Img(:,1:x-1,:);
Gx = k*(xShifted-Img);

grad = sqrt(Gx.^2+Gy.^2+Gz.^2);

end

function [laplacian3D] = find3DLaplacian(Img)

h = zeros(3,3,3);
h(:,:,1) = [0 3 0;3 10 3;0 3 0];
h(:,:,3) = h(:,:,1);
h(:,:,2) = [3 10 3;10 -96 10;3 10 3];

% h(:,:,1) = [0 0 0;0 1 0;0 0 0];
% h(:,:,3) = h(:,:,1);
% h(:,:,2) = [0 1 0;0 -6 0;0 1 0];

k = 1/27;

laplacian3D = k*imfilter(Img,h,'conv','same');

end

