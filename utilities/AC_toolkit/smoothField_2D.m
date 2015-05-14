function [F_hat] = smoothField_2D(F,s )
%SMOOTHFIELD Smooth a vector field by a gaussian and return smooth the unit vector
%field

epsilon = 1e-8;

smoother = fspecial('gaussian',round(7*s),s);

f1 = conv2(F(:,:,1),smoother,'same');
f2 = conv2(F(:,:,2),smoother,'same');

% f_mag = sqrt(f1.^2 + f2.^2) + epsilon;
f_mag = 1;

F_hat(:,:,1) = f1./f_mag;
F_hat(:,:,2) = f2./f_mag;




end

