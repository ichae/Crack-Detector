function [F1,F2,F3,H, enhImg] = medialVectorField3D(Img,smin,smax)
%MEDIALVECTORFIELD3D Obtain the 3D medial flow vector field
%   F1                  : The medial Flow vector field, a 4D array, each
%                         voxel has (x,y,z) components
%   H                   : Hessian matrix at the correct scale, 3x3 matrix
%                         for each voxel linearized.
%   enhImg              : 3D Vessel enhanecd image,  size (row x col x depth)


[r c d]                         = size(Img);
gamma                           = 3;
enhI                            = zeros(r,c,d,(smax-smin+1));
scales                          = smin:smax ;

parfor s = 1 : length(scales)
    sz = 7*scales(s);
    h = fspecial('gaussian',sz, scales(s));
    I_smooth = convn(Img,h,'same');
    [scale_space_hessian(s).hessianMat, enhI(:,:,:,s),eigenVector1(:,:,:,s),eigenVector2(:,:,:,s),eigenVector3(:,:,:,s)] = ...
                                                                                    vesselEnhance3D(I_smooth,scales(s),gamma);
    fprintf('Vesselness for Scale = %d \n', scales(s));
end

% COMMENT:      Compute the vesselness over the scale space

final_enh_I                     = zeros(r,c ,d);


fprintf('Computing the scale space response \n');

for i = 1 : r
    for j = 1 : c
        for k = 1 : d
            
            [final_enh_I(i,j,k),ind] = max(enhI(i,j,k,:));
            finalEV1{i,j,k} = eigenVector1{i,j,k,ind};  % correspond to least eigval
            finalEV2{i,j,k} = eigenVector2{i,j,k,ind};  % correspond to least eigval
            finalEV3{i,j,k} = eigenVector3{i,j,k,ind};  % correspond to least eigval
            %         finalEV2{i,j} = eigenVector2{i,j,ind};  % correspond to greatest eigval
            %         finalEigVal{i,j} = eigenVals{i,j,ind};
            
            linIdx = sub2ind(size(final_enh_I),i,j,k);
            optimum_scale = ind;
            
            hessian_at_optimum_scale = scale_space_hessian(optimum_scale).hessianMat;
            H(:,:,linIdx) = hessian_at_optimum_scale(:,:,linIdx);
            
        end
    end
end

enhImg = final_enh_I;

fprintf ('Stacking up the outputs \n');

for i = 1 : r
    for j = 1 : c
        for k = 1 : d
            v1 = finalEV1{i,j,k};
            F1(i,j,k,1) = v1(1);
            F1(i,j,k,2) = v1(2);
            F1(i,j,k,3) = v1(3);
            
            v2 = finalEV2{i,j,k};
            F2(i,j,k,1) = v2(1);
            F2(i,j,k,2) = v2(2);
            F2(i,j,k,3) = v2(3);
            
            v3 = finalEV3{i,j,k};
            F3(i,j,k,1) = v3(1);
            F3(i,j,k,2) = v3(2);
            F3(i,j,k,3) = v3(3);
            
        end
    end
end



end

