function [F1,enhImg] = obtainVesselTangentField(Img,smin,smax)
%OBTAINVESSELTANGENTFIELD 
% produce the hessian enhanced image: enhI
% return the vessel-flow vector field

[r c] = size(Img);

gamma = 3;
enhI = zeros(r,c,(smax-smin+1));

% scale_space_hessian

for s = smin:smax
    sz = 7*s;
    h = fspecial('gaussian',sz, s);
    I_smooth = conv2(Img,h,'same');
    [hessian, enhI(:,:,s),eigenVector1(:,:,s),eigenVector2(:,:,s),eigenVals(:,:,s)] = vesselEnhance(I_smooth,s,gamma);
    scale_space_hessian.hessianMat(s) = hessian;
end

final_enh_I = zeros(r,c);

for i = 1 : r
    for j = 1 : c
        [final_enh_I(i,j),ind] = max(enhI(i,j,:));
        finalEV1{i,j} = eigenVector1{i,j,ind};  % correspond to least eigval
        finalEV2{i,j} = eigenVector2{i,j,ind};  % correspond to greatest eigval
        finalEigVal{i,j} = eigenVals{i,j,ind};
    end
end

enhImg = final_enh_I;

for i = 1 : r
    for j = 1 : c
        v = finalEV1{i,j};
        F1(i,j,1) = v(1);
        F1(i,j,2) = v(2);
    end
end


% F1 = finalEV1;

end

