function [ G] = createGraph3D_new(Img)
%CREATEGRAPH3D_NEW Summary of this function goes here
%   Detailed explanation goes here


[row col depth] = size(Img);
M = [];
N = [];
W = [];

for i = 2 : row-1
    for j = 2 : col-1
        for k = 2 : depth-1
            current = sub2ind(size(Img),i,j,k);
            % look for 6 neighbours initially
            n1 = sub2ind(size(Img),i+1,j,k);
            n2 = sub2ind(size(Img),i-1,j,k);
            n3 = sub2ind(size(Img),i,j+1,k);
            n4 = sub2ind(size(Img),i,j-1,k);
            n5 = sub2ind(size(Img),i,j,k+1);
            n6 = sub2ind(size(Img),i,j,k-1);
            
            n7 = sub2ind(size(Img),i+1,j+1,k);
            n8 = sub2ind(size(Img),i+1,j-1,k);
            n9 = sub2ind(size(Img),i-1,j+1,k);
            n10 = sub2ind(size(Img),i-1,j-1,k);
            
            n11 = sub2ind(size(Img),i+1,j+1,k+1);
            n12 = sub2ind(size(Img),i+1,j+1,k-1);
            
            n13 = sub2ind(size(Img),i+1,j-1,k+1);
            n14 = sub2ind(size(Img),i+1,j-1,k-1);
            
            n15 = sub2ind(size(Img),i-1,j+1,k+1);
            n16 = sub2ind(size(Img),i-1,j+1,k-1);
            
            n17 = sub2ind(size(Img),i-1,j-1,k+1);
            n18 = sub2ind(size(Img),i-1,j-1,k-1);
            
            n19 = sub2ind(size(Img),i+1,j,k+1);
            n20 = sub2ind(size(Img),i+1,j,k-1);
            
            n21 = sub2ind(size(Img),i-1,j,k+1);
            n22 = sub2ind(size(Img),i-1,j,k-1);
            
            n23 = sub2ind(size(Img),i,j-1,k+1);
            n24 = sub2ind(size(Img),i,j-1,k-1);
            
            n25 = sub2ind(size(Img),i,j+1,k+1);
            n26 = sub2ind(size(Img),i,j+1,k-1);
            
            
            
            M = [M ;current]; N = [N;n1]; W = [W;edgeWeight(i,j,k,i+1,j,k,Img)];
            M = [M;n1]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j,k,Img)];
            
            M = [M ;current]; N = [N;n2]; W = [W;edgeWeight(i,j,k,i-1,j,k,Img)];
            M = [M;n2]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j,k,Img)];
            
            M = [M ;current]; N = [N;n3]; W = [W;edgeWeight(i,j,k,i,j+1,k,Img)];
            M = [M;n3]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j+1,k,Img)];
            
            M = [M ;current]; N = [N;n4]; W = [W;edgeWeight(i,j,k,i,j-1,k,Img)];
            M = [M;n4]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j-1,k,Img)];
            
            M = [M ;current]; N = [N;n5]; W = [W;edgeWeight(i,j,k,i,j,k+1,Img)];
            M = [M;n5]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j,k+1,Img)];
            
            M = [M ;current]; N = [N;n6]; W = [W;edgeWeight(i,j,k,i,j,k-1,Img)];
            M = [M;n6]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j,k-1,Img)];
            
            M = [M ;current]; N = [N;n7]; W = [W;edgeWeight(i,j,k,i+1,j+1,k,Img)];
            M = [M;n7]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j+1,k,Img)];
            
            M = [M ;current]; N = [N;n8]; W = [W;edgeWeight(i,j,k,i+1,j-1,k,Img)];
            M = [M;n8]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j-1,k,Img)];
            
            M = [M ;current]; N = [N;n9]; W = [W;edgeWeight(i,j,k,i-1,j+1,k,Img)];
            M = [M;n9]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j+1,k,Img)];
            
            M = [M ;current]; N = [N;n10]; W = [W;edgeWeight(i,j,k,i-1,j-1,k,Img)];
            M = [M;n10]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j-1,k,Img)];
            
            M = [M ;current]; N = [N;n11]; W = [W;edgeWeight(i,j,k,i+1,j+1,k+1,Img)];
            M = [M;n11]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j+1,k+1,Img)];
            
            M = [M ;current]; N = [N;n12]; W = [W;edgeWeight(i,j,k,i+1,j+1,k-1,Img)];
            M = [M;n12]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j+1,k-1,Img)];
            
            M = [M ;current]; N = [N;n13]; W = [W;edgeWeight(i,j,k,i+1,j-1,k+1,Img)];
            M = [M;n13]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j-1,k+1,Img)];
            
            M = [M ;current]; N = [N;n14]; W = [W;edgeWeight(i,j,k,i+1,j-1,k-1,Img)];
            M = [M;n14]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j-1,k-1,Img)];
            
            M = [M ;current]; N = [N;n15]; W = [W;edgeWeight(i,j,k,i-1,j+1,k+1,Img)];
            M = [M;n15]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j+1,k+1,Img)];
            
            M = [M ;current]; N = [N;n16]; W = [W;edgeWeight(i,j,k,i-1,j+1,k-1,Img)];
            M = [M;n16]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j+1,k-1,Img)];
            
            M = [M ;current]; N = [N;n17]; W = [W;edgeWeight(i,j,k,i-1,j-1,k+1,Img)];
            M = [M;n17]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j-1,k+1,Img)];
            
            M = [M ;current]; N = [N;n18]; W = [W;edgeWeight(i,j,k,i-1,j-1,k-1,Img)];
            M = [M;n18]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j-1,k-1,Img)];
            
            M = [M ;current]; N = [N;n19]; W = [W;edgeWeight(i,j,k,i+1,j,k+1,Img)];
            M = [M;n19]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j,k+1,Img)];
            
            M = [M ;current]; N = [N;n20]; W = [W;edgeWeight(i,j,k,i+1,j,k-1,Img)];
            M = [M;n20]; N = [N;current];  W = [W;edgeWeight(i,j,k,i+1,j,k-1,Img)];
            
            M = [M ;current]; N = [N;n21]; W = [W;edgeWeight(i,j,k,i-1,j,k+1,Img)];
            M = [M;n21]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j,k+1,Img)];
            
            M = [M ;current]; N = [N;n22]; W = [W;edgeWeight(i,j,k,i-1,j,k-1,Img)];
            M = [M;n22]; N = [N;current];  W = [W;edgeWeight(i,j,k,i-1,j,k-1,Img)];
            
            M = [M ;current]; N = [N;n23]; W = [W;edgeWeight(i,j,k,i,j-1,k+1,Img)];
            M = [M;n23]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j-1,k+1,Img)];
            
            M = [M ;current]; N = [N;n24]; W = [W;edgeWeight(i,j,k,i,j-1,k-1,Img)];
            M = [M;n24]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j-1,k-1,Img)];
            
            M = [M ;current]; N = [N;n25]; W = [W;edgeWeight(i,j,k,i,j+1,k+1,Img)];
            M = [M;n25]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j+1,k+1,Img)];
            
            M = [M ;current]; N = [N;n26]; W = [W;edgeWeight(i,j,k,i,j+1,k-1,Img)];
            M = [M;n26]; N = [N;current];  W = [W;edgeWeight(i,j,k,i,j+1,k-1,Img)];
        end
    end
end

r = row*col*depth;
G = sparse(M,N,W,r,r);

end


function wt = edgeWeight(xp,yp,zp,xq,yq,zq,Img)

lambda = .001;
wt = 1/(lambda+abs(Img(xp,yp,zp) + Img(xq,yq,zq)));

end

