function outimage=areaopenfn2(x,minarea)
%  Area open filter

[height, width]=size(x);
x = x*255;
figure; imshow(x,[]);      
L = zeros(height, width, 255);  % stack of level sets
for l=1:255,
   levelset=(x>=l);
   levelset=bwareaopen(levelset,minarea);
   imshow(~levelset);drawnow;
   L(:,:,l)=levelset;
end
outimage=sum(L,3);
outimage = outimage/255;




