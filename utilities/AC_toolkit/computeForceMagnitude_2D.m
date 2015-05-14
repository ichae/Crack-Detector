function [M] = computeForceMagnitude_2D(F,I,enhI,tau)
%COMPUTEFORCEMAGNITUDE_Compute the force magnitude based on the directional
%derivative of the normal vector.

sigma = 2;
h = fspecial('gaussian',7*sigma,sigma);
I = conv2(I,h,'same');
[gy,gx] = gradient(I);


grad_mat(:,:,1) = gy;
grad_mat(:,:,2) = gx;

directional_gradient = dot(F,grad_mat,3);

M = 0.5*(enhI + exp(-abs(directional_gradient./(I + 0.001))/tau));

M = (M - min(M(:)))/(max(M(:))-min(M(:)));


end

