function [MIP] = getMIP(Img)
%SHOWMIP display the Maximum intensity Projection

[row col depth] = size(Img);

MIP = zeros(row,col);

for i = 1 : row
    for j = 1 : col
        MIP(i,j) = max(Img(i,j,:));
    end
end

end

