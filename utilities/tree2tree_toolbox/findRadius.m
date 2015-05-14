function [finalRad] = findRadius(Img, x , y, z ,Rmin, Rmax)
%FINDRADIUS Find the actual value of the radii at a neuronal point


% Rmin = 1 ;
% Rmax = 4 ;

interval = 0.5 ;


index = 1;
pValue = findIntensity(x,y,z, Rmin, Img);

for R = Rmin : interval : Rmax
    
    value(index) = findIntensity(x,y,z, R , Img);
    diff(index) = abs(value(index) - pValue);
    pValue = value(index);
    currentRadii(index) = R ;
    index = index+1;
end

% disp(value');
% pause;

[val ind] = max(diff); % find maximum rate of change
finalRad = currentRadii(ind);

end



function [intensity] = findIntensity(x_c, y_c, z_c, R, Img)

angleInterval = 20 ;
intensity = 0;
[row col depth] = size(Img);

for theta = -pi/2 : pi/angleInterval : pi/2
   for phi = 0 :   pi/angleInterval  : pi
       x = x_c+ round(R*cos(theta)*cos(phi)); 
       y = y_c+ round(R*sin(theta)*cos(phi));
       z = z_c+ round(R*sin(phi));
       
       if (x >= 1 && y >= 1 && z >= 1 && x <= row ...
               && y <= col && z <= depth)
           intensity = intensity + Img(x,y,z);
       end
   end
end

% intensity = intensity / (4*pi*R*R) ;
end