function [ hessian , enhI, evCell1, evCell2, evCell3 ] = vesselEnhance3D(I,s,gamma)
%VESSELENHANCE3D Compute the 3D vesselness for the scale s
% hessian:  (3*3)x(row x col x depth) matrix, i.e a 3*3 matrix for each voxel
% enhI   :  enhanced Image  , size = row x col x depth
% evCell1:  eigencevtor for each voxel, corresponding to lowest eigenvalue. size: cell( linear(row x col x depth)). Each cell has 3 components x,y,z  

[nrow ncol nslice]               = size(I);

[dIdX , dIdY,  dIdZ]             = gradient(I);
[d2IdX2, d2IdXdY, d2IdXdZ ]      = gradient(dIdX);
[~, d2IdY2, d2IdYdZ]             = gradient(dIdY);
[~, ~ , d2IdZ2]                  = gradient(dIdZ);

% Compute the Hessian
hessian                          = zeros(3,3,nrow*ncol*nslice) ;

hessian(1,1,:)                   = d2IdX2(:);
hessian(2,2,:)                   = d2IdY2(:);
hessian(3,3,:)                   = d2IdZ2(:);

hessian(1,2,:)                   = d2IdXdY(:);
hessian(2,1,:)                   = hessian(1,2,:);

hessian(1,3,:)                   = d2IdXdZ(:);
hessian(3,1,:)                   = hessian(1,3,:);

hessian(2,3,:)                   = d2IdYdZ(:);
hessian(3,2,:)                   = hessian(2,3,:);

enhI                             = zeros(length(I(:)),1) ;
alpha                            = 0.5;
beta                             = 0.5;
c                                = 0.5;
epsilon                          = 0.0;   

len = length(I(:));

for i = 1 : len
    tempHess = hessian(:,:,i);
    [V,e] = eig(tempHess);
    D                           = diag(e);         % store the signed eigenvalues in D 
    [~,sortindex]               = sort(abs(D));    % sorted in ascending order 
    evCell1{i}                  = [0;0;0];         % eigenvector along medial axis
    evCell2{i}                  = [0;0;0];         % eigenvector along medial axis
    evCell3{i}                  = [0;0;0];         % eigenvector along medial axis
    
    if (D(sortindex(2)) < 0 && D(sortindex(3)) < 0)% if the larger eigenvalue is negative, for bright object

        %COMMENT:  |l1| < = |l2|  <= |l3|
        
        l1                      = D(sortindex(1));
        l2                      = D(sortindex(2));
        l3                      = D(sortindex(3));
        
        R_A = abs(l2)/(abs(l3) + epsilon);
        R_B = abs(l1)/(epsilon + sqrt(abs(l2)*abs(l3)));
        S = sqrt(l1^2+l2^2+l3^2);
        
        vesselness = (1 - exp(-(R_A^2)/(2*alpha^2)))*(exp(-(R_B^2)/(2*beta^2)))*(1-exp(-(S^2)/(2*c^2)));
        evCell1{i} = V(:,sortindex(1));
        evCell2{i} = V(:,sortindex(2));
        evCell3{i} = V(:,sortindex(3));
        enhI(i) = vesselness ;
        
    end
end

enhI = enhI*s^gamma;
enhI = reshape(enhI,size(I));
evCell1 = reshape(evCell1,size(I));
evCell2 = reshape(evCell2,size(I));
evCell3 = reshape(evCell3,size(I));


end

