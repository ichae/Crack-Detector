function [attraction_force] = attractionForce(binaryImg,gamma,epsilon)
%ATTRACTIONFORCE Calculate the force of attraction using the VFC
% binaryImg: phi >=0;
% gamma    : VFC kernel parameter


Dirac                           = @(x,e) ((1/2/e)*(1+cos(pi*x/e)).*((x<=e) & (x>=-e)));
Heaviside                       = @(y,e) (0.5*(1+(2/pi)*atan(y/e)));

CC = bwconncomp(binaryImg);
num_cc = CC.NumObjects;

if num_cc <= 1
    attraction_force = zeros(size(binaryImg));
else
    
    rp = regionprops(CC, 'Area', 'PixelIdxList');
    [area,ind] = max([rp.Area]);
    rp = rp(ind);
    
    largest_bwI = false(size(binaryImg));
    largest_bwI(rp(1).PixelIdxList) = true;
    
    
    % Compute the SDF for only the largest CC
    phi_single          = double(bwdist(~largest_bwI) - bwdist(largest_bwI));
    
    delta_phi_largest   = Dirac(phi_single,epsilon);
    delta_phi_largest   = (delta_phi_largest - min(delta_phi_largest(:)))/(max(delta_phi_largest(:))-min(delta_phi_largest(:)));

    
    % Compute the SDF for the remaining components
    
    other_bwI = binaryImg - largest_bwI;
    surface_map =  double(bwdist(~other_bwI) - bwdist(other_bwI));
    
    % Compute the attraction field using VFC

    rad = max(round(3*gamma),1);
    Ker = AM_VFK(2,rad,'gaussian',gamma);   % Gaussian kernel of radii = rad
    
    vfkX = Ker(:,:,1);
    vfkY = Ker(:,:,2);
    
    Fx = conv2(delta_phi_largest,vfkX,'same').*(~largest_bwI);
    Fy = conv2(delta_phi_largest,vfkY,'same').*(~largest_bwI);
    
    
    % Compute the attraction force
    
    [gx,gy]             = gradient(surface_map);
    mag                 = sqrt(gx.^2+gy.^2 + 1e-10);
    kappa               = -divergence(gx./mag,gy./mag);
        
    delta_phi           = Dirac(surface_map,10);
    outward_normal(:,:,1)      = -gx.*delta_phi;
    outward_normal(:,:,2)      = -gy.*delta_phi;     % negative sign because LSF is positive inside, hence concave
    attraction_field(:,:,1)    = Fx;
    attraction_field(:,:,2)    = Fy;
    
    Mass = 1;
    attraction_force = Mass*dot(outward_normal,attraction_field,3);
    attraction_force = attraction_force.*kappa;

    if max(attraction_force(:) ~=0)
        disp('attraction');
        attraction_force = (attraction_force - min(attraction_force(:)))/(max(attraction_force(:))-min(attraction_force(:)));
    end
end

end

