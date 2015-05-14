function enhI = darkVesselEnhanceFrangi(I,scales,gamma)
% Enhance Bright vessels

[r,c] = size(I);
n = length(scales);
enh_stack = zeros(r,c,n);
for ii = 1 : length(scales)
    sigma = scales(ii);
    sz = 6*sigma;
    h = fspecial('gaussian',sz, sigma);
    I_smooth = conv2(I,h,'same');
    frangi_resp = sigma^gamma*vesselEnhanceDark(I_smooth);
    enh_stack(:,:,ii)=frangi_resp;
end
enhI = max(enh_stack,[],3);
end



function [enhI] = vesselEnhanceDark(I)
%VESSELENHANCE Summary of this function goes here
%   Detailed explanation goes here

[row, col]= size(I);

[dIdX, dIdY] = gradient(I);
[d2IdX2, d2IdXdY] = gradient(dIdX);
[d2IdYdX, d2IdY2] = gradient(dIdY);

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

prevVesselness = 0 ;

len = length(I(:));


for i = 1 : len
    tempHess = hessian(:,:,i);
    e = eig(tempHess);
    [eigVal,sortindex]=sort(abs(e));
    e1 = eigVal(1);
    e2 = eigVal(2);
    S = sqrt(e1^2 + e2^2);
    if (e(sortindex(2)) > 0) % if the smaller eigenvalue is negative
        vesselness = findVesselness(e(sortindex(1)),e(sortindex(2)),S,c,beta);
        enhI(i) = vesselness ;
    end
end
enhI = reshape(enhI,size(I));
enhI = (enhI-min(enhI(:)))/(max(enhI(:))-min(enhI(:)));
% enhI = enhI./max(enhI(:));

end

function [v] = findVesselness(e1,e2,s,c,beta)
%Find the vesselness cost function

if (e2 < 0)
    v = 0 ;
    disp('eshechi');
else
   Rb = abs(e1/(e2+0.0001)); % blobness measure
   v = exp(-Rb^2/(2*beta^2))*(1 - exp(-s^2/(2*c^2)));
end

end
