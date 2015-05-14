function [ G ] = createGraph3D(Img )
%CREATEGRAPH3D Summary of this function goes here


[row col depth] = size(Img);

% G = zeros(row*col*depth);
M = []; N = []; W = [];
index = 1 ;
for i = 2 : row-1
    for j = 2 : col-1
        for k = 2 : depth-1
            
            M(index,1) = sub2ind(size(Img),i,j,k); N(index,1) = sub2ind(size(Img),i,j+1,k) ; W(index,1) = edgeWeight(i,j,k,i,j+1,k,Img);
            N(index+1,1) = sub2ind(size(Img),i,j,k); M(index+1,1) = sub2ind(size(Img),i,j+1,k) ; W(index+1,1) = edgeWeight(i,j,k,i,j+1,k,Img);
            
            M(index+2,1) = sub2ind(size(Img),i,j,k); N(index+2,1) = sub2ind(size(Img),i+1,j,k) ; W(index+2,1) = edgeWeight(i,j,k,i+1,j,k,Img);
            N(index+3,1) = sub2ind(size(Img),i,j,k); M(index+3,1) = sub2ind(size(Img),i+1,j,k) ; W(index+3,1) = edgeWeight(i,j,k,i+1,j,k,Img);
            
            M(index+4,1) = sub2ind(size(Img),i,j,k); N(index+4,1) = sub2ind(size(Img),i+1,j+1,k) ; W(index+4,1) = edgeWeight(i,j,k,i+1,j+1,k,Img);
            N(index+5,1) = sub2ind(size(Img),i,j,k); M(index+5,1) = sub2ind(size(Img),i+1,j+1,k) ; W(index+5,1) = edgeWeight(i,j,k,i+1,j+1,k,Img);
            
            M(index+6,1) = sub2ind(size(Img),i,j,k); N(index+6,1) = sub2ind(size(Img),i,j,k+1) ; W(index+6,1) = edgeWeight(i,j,k,i,j,k+1,Img);
            M(index+7,1) = sub2ind(size(Img),i,j,k); N(index+7,1) = sub2ind(size(Img),i,j,k+1) ; W(index+7,1) = edgeWeight(i,j,k,i,j,k+1,Img);

            index = index+8;    
%             G(sub2ind(size(Img),i,j,k),sub2ind(size(Img),i,j+1,k)) = edgeWeight(i,j,k,i,j+1,k,Img);
%             G(sub2ind(size(Img),i,j,k),sub2ind(size(Img),i+1,j,k)) = edgeWeight(i,j,k,i+1,j,k,Img);
%             G(sub2ind(size(Img),i,j,k),sub2ind(size(Img),i+1,j+1,k)) = edgeWeight(i,j,k,i+1,j+1,k,Img);
%             G(sub2ind(size(Img),i,j,k),sub2ind(size(Img),i,j,k+1)) = edgeWeight(i,j,k,i,j,k+1,Img);
        end
    end
end
% M
% r = max(max(M(:)),max(N(:)));
r = row*col*depth;
G = sparse(M,N,W,r,r);


end



function wt = edgeWeight(xp,yp,zp,xq,yq,zq,Img)

lambda = 1;
wt = 1/(lambda+abs(Img(xp,yp,zp) + Img(xq,yq,zq)));

end



