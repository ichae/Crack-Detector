clear all;
clc;

% folder = 

num_files = 239; % change here

figure ; hold on;

filename = 'animated.gif';
for ind = 1 : num_files
    fname = ['a' num2str(ind) '.BMP'];
    Img = imread(fname);
    imagesc(Img);colormap('gray');
    
     frame = getframe(1);
     im = frame2im(frame);
     [imind,cm] = rgb2ind(im,256);
     if ind == 1;
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
     else
         imwrite(imind,cm,filename,'gif','WriteMode','append');
     end
    
end

hold off
% for ind = 1 : num_files
%    
%     fname = num2str(ind);
%     fullfname = [fname '.bmp'];
%     I = imread(fullfname);
%     newname = [fname '.gif'];
%     imwrite(I,newname);
%     fname = [];
%     fullfname = [];
%     newname = [];
% end