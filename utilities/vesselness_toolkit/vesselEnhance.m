function [enhI, evCell1,evCell2, evalCell] = vesselEnhance(I,s,gamma)
%VESSELENHANCE Summary of this function goes here
%   Detailed explanation goes here

[row col]= size(I);

[dIdX dIdY] = gradient(I);
[d2IdX2 d2IdXdY] = gradient(dIdX);
[d2IdYdX d2IdY2] = gradient(dIdY);

disp('vessel Enhance');

hessian = zeros(2,2,row*col) ;

hessian(1,1,:) = d2IdX2(:);
hessian(2,2,:) = d2IdY2(:);
hessian(1,2,:) = d2IdXdY(:);
hessian(2,1,:) = hessian(1,2,:);

tempHess = zeros(2,2);

count = 0;
enhI = zeros(length(I(:)),1) ;

% beta and c are two controlling parameters
beta = 0.5;
c = 0.3 ;
Smin = 0.2 ;
Smax = 0.8;

prevVesselness = 0 ;

len = length(I(:));


for i = 1 : len
    tempHess = hessian(:,:,i);
    [V,e] = eig(tempHess);
    [eigVal,sortindex]=sort(abs(diag(e)));
    e1 = eigVal(1);
    e2 = eigVal(2);
    
    evCell1{i} = [0;0];
    evCell2{i} = [0;0];
    evalCell{i} = [e1;e2];

    S = sqrt(e1^2 + e2^2);
    if (e(sortindex(2))<=0) % if the smaller eigenvalue is negative, for bright object
        vesselness = findVesselness(e(sortindex(1)),e(sortindex(2)),S,c,beta);
        evCell1{i} = V(:,sortindex(1));
        evCell2{i} = V(:,sortindex(2));
        enhI(i) = vesselness ;
    end
end
enhI = enhI*s^gamma;
enhI = reshape(enhI,size(I));

evCell1 = reshape(evCell1,size(I));
evCell2 = reshape(evCell2,size(I));
evalCell = reshape(evalCell,size(I));
end





