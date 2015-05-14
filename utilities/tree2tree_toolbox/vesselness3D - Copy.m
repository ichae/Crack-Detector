
function vslI = vesselness3D(Img)
% Computes the vesselness measure of the 3D image 'Img'
% 'Img' is smoothed by a Gaussian of size 's' before calling this function

vslI=zeros(size(Img(:)));

z_pshifted = circshift(Img,[0 0 1]);
z_nshifted = circshift(Img,[0 0 -1]);
x_pshifted = circshift(Img,[0 1 0]);
x_nshifted = circshift(Img,[0 -1 0]);
y_pshifted = circshift(Img,[1 0 0]);
y_nshifted = circshift(Img,[-1 0 0]);

dIx=x_pshifted-x_nshifted;
dIy=y_pshifted-y_nshifted;
dIz=z_pshifted-z_nshifted;

dIx_pshifted=circshift(dIx,[0 1 0]);
dIx_nshifted=circshift(dIx,[0 -1 0]);
dIy_pshifted=circshift(dIy,[1 0 0]);
dIy_nshifted=circshift(dIy,[-1 0 0]);
dIz_pshifted=circshift(dIz,[0 0 1]);
dIz_nshifted=circshift(dIz,[0 0 -1]);

d2Idx2=dIx_pshifted-dIx_nshifted;
d2Idx2=d2Idx2(:);
d2Idy2=dIy_pshifted-dIy_nshifted;
d2Idy2=d2Idy2(:);
d2Idz2=dIz_pshifted-dIz_nshifted;
d2Idz2=d2Idz2(:);

dIx_pyshifted=circshift(dIx,[1 0 0]);
dIx_nyshifted=circshift(dIx,[-1 0 0]);
d2Idxdy=dIx_pyshifted-dIx_nyshifted;
d2Idxdy=d2Idxdy(:);

dIy_pzshifted=circshift(dIy,[0 0 1]);
dIy_nzshifted=circshift(dIy,[0 0 -1]);
d2Idydz=dIy_pzshifted-dIy_nzshifted;
d2Idydz=d2Idydz(:);

dIz_pxshifted=circshift(dIz,[0 1 0]);
dIz_nxshifted=circshift(dIz,[0 -1 0]);
d2Idzdx=dIz_pxshifted-dIz_nxshifted;
d2Idzdx=d2Idzdx(:);

npxls=length(Img(:));

HessianStack = zeros(3,3,npxls);

HessianStack(1,1,:) = d2Idx2;
HessianStack(2,2,:) = d2Idy2;
HessianStack(3,3,:) = d2Idz2;
HessianStack(1,2,:) = d2Idxdy;
HessianStack(2,1,:) = d2Idxdy;
HessianStack(1,3,:) = d2Idzdx;
HessianStack(3,1,:) = d2Idzdx;
HessianStack(2,3,:) = d2Idydz;
HessianStack(3,2,:) = d2Idydz;

        alpha = 0.5;
        beta = 0.5;
        c = 0.5;
for I = 1 : npxls,
%     I
    % take the eigenvalue for each pixel
    testmat = HessianStack(:,:,I);
    [V,D] = eig(testmat);
    D = diag(D);
    [eigVal,sortindex] = sort(abs(D));
    
    
    if D(sortindex(2))<0 && D(sortindex(3))<0
        l1 = D(sortindex(1));
        l2 = D(sortindex(2));
        l3 = D(sortindex(3));
        
        R_A = abs(l2)/abs(l3);
        R_B = abs(l1)/sqrt(abs(l2)*abs(l3));
        S = sqrt(l1^2+l2^2+l3^2);
        
        vslI(I) = (1 - exp(-(R_A^2)/(2*alpha^2)))*(exp(-(R_B^2)/(2*beta^2)))*(1-exp(-(S^2)/(2*c^2)));
        if isnan(vslI(I))
            fprintf('\n NAN');
            disp([R_A R_B S]);
        end
%         vslI(I)=(eigVal(2)-eigVal(1))^2/(abs(eigVal(1))*abs(eigVal(3)-eigVal(2))+0.001);
        
    end
    
%     vslI(I)=(eigVal(2)-eigVal(1))/(eigVal(1)*abs(eigVal(2)-eigVal(3)));
        
%     eigVec{I}=V(:,sortindex);
    
end

% vslI=reshape(vslI,size(Img));


