function [] = compareVisually(I,resp,frangi_resp,thresh_red,spur_amt,area,mag)
%COMPAREVISUALLY Compare the results of Frangi and LDE
% Plot enhancement result
figure(4);
subtightplot(1,2,1,0.01,0.01,0.01); imshow(resp,[],'InitialMagnification',mag); title('Our');
subtightplot(1,2,2,0.01,0.01,0.01); imshow(frangi_resp,[],'InitialMagnification',mag); title('Frangi');


% Plot the binarized image
bwI = resp > thresh_red*graythresh(resp);
bwI = bwareaopen(bwI,area);
se = strel('disk',1);
% bwI = imclose(bwI,se);

bwI_f = frangi_resp > thresh_red*graythresh(frangi_resp);
bwI_f = bwareaopen(bwI_f,area);
se = strel('disk',1);
bwI_f = imclose(bwI_f,se);

figure(5); 
subtightplot(1,2,1,0.01,0.01,0.01);imshow(bwI,'InitialMagnification',mag); title('Ours');
subtightplot(1,2,2,0.01,0.01,0.01);imshow(bwI_f,'InitialMagnification',mag); title('Frangi');

% Plot the skeleton
skel = bwmorph(bwI,'skel',Inf);
skel = bwmorph(skel,'spur',spur_amt);
[R,C] = find(skel);
figure(6); 
subtightplot(1,2,1,0.01,0.01,0.01);imshow(I,[],'InitialMagnification',mag); title('Ours');hold on;
plot(C,R,'.g','MarkerSize',5); 
hold off;

skel_f = bwmorph(bwI_f,'skel',Inf);
skel_f = bwmorph(skel_f,'spur',spur_amt);
[R,C] = find(skel_f);
figure(6); 
subtightplot(1,2,2,0.01,0.01,0.01);imshow(I,[],'InitialMagnification',mag);title('frangi'); hold on;
plot(C,R,'.c','MarkerSize',5); 
hold off;



end

