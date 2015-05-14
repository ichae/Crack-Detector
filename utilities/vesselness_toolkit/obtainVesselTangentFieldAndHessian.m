function [F1,H, enhImg] = obtainVesselTangentFieldAndHessian(Img,smin,smax)
%OBTAINVESSELTANGENTFIELD 
% produce the hessian enhanced image: enhI
% return the vessel-flow vector field
% return the optimum scale hessian for each pixel. The pixels are indexed
% linearly.(Matlab follows column major order)

[r c] = size(Img);

gamma = 3;
enhI = zeros(r,c,(smax-smin+1));

% scale_space_hessian
scales = smin:smax;
for s = 1:length(scales)
    sz = 7*scales(s);
    h = fspecial('gaussian',sz, scales(s));
    I_smooth = conv2(Img,h,'same');
    [hessian, enhI(:,:,s),eigenVector1(:,:,s),eigenVector2(:,:,s),eigenVals(:,:,s)] = vesselEnhance2(I_smooth,scales(s),gamma);
    scale_space_hessian(s).hessianMat = hessian;
end

final_enh_I = zeros(r,c);

for i = 1 : r
    for j = 1 : c
        [final_enh_I(i,j),ind] = max(enhI(i,j,:));
        finalEV1{i,j} = eigenVector1{i,j,ind};  % correspond to least eigval
        finalEV2{i,j} = eigenVector2{i,j,ind};  % correspond to greatest eigval
        finalEigVal{i,j} = eigenVals{i,j,ind};
        
        linIdx = sub2ind(size(final_enh_I),i,j);
        optimum_scale = ind;
        
        hessian_at_optimum_scale = scale_space_hessian(optimum_scale).hessianMat;
        H(:,:,linIdx) = hessian_at_optimum_scale(:,:,linIdx);
        
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

