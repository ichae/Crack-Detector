function [F_total ] = attractionForce3D( binaryImg,gamma,capture_range )
%ATTRACTIONFORCE3D 3D version of the attraction force

Dirac                           = @(x,e) ((1/2/e)*(1+cos(pi*x/e)).*((x<=e) & (x>=-e)));
Heaviside                       = @(y,e) (0.5*(1+(2/pi)*atan(y/e)));
num_masters     = 2;


CC = bwconncomp(binaryImg);
num_cc = CC.NumObjects;
F_total = zeros(size(binaryImg));
if num_cc <= 1
    attraction_force = zeros(size(binaryImg));
else
    
    regions = regionprops(CC, 'Area', 'PixelIdxList');
    %     [area,ind] = max([rp.Area]);
    [~,sortInd] = sort([regions.Area],'descend');
    
    for ii = 1 : num_masters
        rp = regions(sortInd(ii));
        
        largest_bwI = false(size(binaryImg));
        largest_bwI(rp(1).PixelIdxList) = true;
        
        % Compute the SDF for only the largest CC
        phi_single          = double(bwdist(~largest_bwI) - bwdist(largest_bwI));
        delta_phi_largest   = Dirac(phi_single,capture_range);
        delta_phi_largest   = (delta_phi_largest - min(delta_phi_largest(:)))/(max(delta_phi_largest(:))-min(delta_phi_largest(:)));
        
        % Compute the SDF for the remaining components
        
        other_bwI = binaryImg - largest_bwI;
        phi_other =  double(bwdist(~other_bwI) - bwdist(other_bwI));
        
        % Compute the attraction field using VFC
        
        rad = max(round(15*capture_range),1);
        Ker = AM_VFK(3,rad,'power',gamma);   % Gaussian kernel of radii = rad
        
        vfkX = Ker(:,:,1);
        vfkY = Ker(:,:,2);
        vfkZ = Ker(:,:,3);
        
        Fx = convn(delta_phi_largest,vfkX,'same').*(~largest_bwI);
        Fy = convn(delta_phi_largest,vfkY,'same').*(~largest_bwI);
        Fz = convn(delta_phi_largest,vfkZ,'same').*(~largest_bwI);
        
        % Compute the attraction force
        
        [gx,gy,gz]              = gradient(phi_other);
        mag                     = sqrt(gx.^2+gy.^2+gz.^2 + 1e-10);
        kappa                   = -divergence(gx./mag,gy./mag,gz./mag);
        
        delta_phi               = Dirac(phi_other,capture_range);
        
        
        nx      = -gx.*delta_phi;
        ny      = -gy.*delta_phi;     % negative sign because LSF is positive inside, hence concave
        nz      = -gz.*delta_phi;
        e_x    = Fx;
        e_y    = Fy;
        e_z    = Fz;
        
        
        Mass = 1;
        attraction_force = Mass*(nx.*e_x+ny.*e_y+nz.*e_z);
        attraction_force = attraction_force.*kappa;
        if max(attraction_force(:) ~=0)
            disp('attraction');
            attraction_force = (attraction_force - min(attraction_force(:)))/(max(attraction_force(:))-min(attraction_force(:)));
        end
        
    end
    F_total = F_total+attraction_force;
end

end

