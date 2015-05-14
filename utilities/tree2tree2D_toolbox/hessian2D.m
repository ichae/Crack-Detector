
function vslI=hessian2D(Img)

Img=double(Img);

[dIdx,dIdy]=gradient(Img);
[d2Idx2,d2Idxdy]=gradient(dIdx);
[d2Idy2,d2Idydx]=gradient(dIdy);
npxls=length(Img(:));

sum(sum(d2Idxdy-d2Idydx))

Hessian=zeros(2,2,npxls);
Hessian(1,1,:)=d2Idx2(:);
Hessian(2,2,:)=d2Idy2(:);
Hessian(1,2,:)=d2Idxdy(:);
Hessian(2,1,:)=Hessian(1,2,:);
vslI=zeros(npxls,1);

for I=1:npxls,
%     I
    testmat=Hessian(:,:,I);
    [V,D] = eig(testmat);
    D=diag(D);
    [eigVal,sortindex]=sort(abs(D));
    
    if D(sortindex(2))<=0,
        
        vslI(I)=(eigVal(2)-eigVal(1))/(abs(eigVal(1))+0.001);
        
    end
    
end

vslI=reshape(vslI,size(Img));