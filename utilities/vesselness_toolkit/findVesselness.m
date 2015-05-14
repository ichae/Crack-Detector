function [v] = findVesselness(e1,e2,s,c,beta)
%Find the vesselness cost function

if (e2 > 0)
    v = 0 ;
    disp('eshechi');
else
   Rb = abs(e1/(e2+0.0001)); % blobness measure
   v = exp(-Rb^2/(2*beta^2))*(1 - exp(-s^2/(2*c^2)));
end
% disp('hello 2')


% v = (e2-e1)/(abs(e1)+.001);


end

