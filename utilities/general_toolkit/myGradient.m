function [Gx,Gy] = myGradient(Img)
%MYGRADIENT find the gradient. Use this instead of matlab's 'gradient'
% Gx is the gradient in X- direction: i.e the direction of rows (top to bottom)
[row col] = size(Img);
xShifted = zeros(size(Img));
yShifted = zeros(size(Img));

xShifted(2:row,:) = Img(1:row-1,:);
yShifted(:,2:col) = Img(:,1:col-1);

Gx = zeros(size(Img));
Gy = zeros(size(Img));

Gx = Img - xShifted ;
Gy = Img - yShifted ;

end

