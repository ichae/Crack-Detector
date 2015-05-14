function [F_total] = attractionForceCurvatureMap(binaryImg,gamma,epsilon)
%ATTRACTIONFORCE Calculate the force of attraction using the VFC
% binaryImg: phi >=0;
% gamma    : VFC kernel parameter

use_convex_hull = 1;
num_masters     = 2;
Dirac                           = @(x,e) ((1/2/e)*(1+cos(pi*x/e)).*((x<=e) & (x>=-e)));
Heaviside                       = @(y,e) (0.5*(1+(2/pi)*atan(y/e)));
% epsilon                         = 2;

CC = bwconncomp(binaryImg);
num_cc = CC.NumObjects;
F_total = zeros(size(binaryImg));

if num_cc <= 1
    attraction_force = zeros(size(binaryImg));
else
    
    regions = regionprops(CC, 'Area', 'PixelIdxList');
    %     [~,ind] = max([rp.Area]);
    [~,sortList] = sort([regions.Area],'descend');
    
    for ii = 1 : num_masters
        rp = regions(sortList(ii));
        
        largest_bwI = false(size(binaryImg));
        largest_bwI(rp(1).PixelIdxList) = true;
        
        
        % Compute the SDF for only the largest CC
        
        phi_single          = double(bwdist(~largest_bwI) - bwdist(largest_bwI));
        CH_largest = bwconvhull(largest_bwI);
        delta_phi_largest   = Dirac(phi_single,epsilon);
%         delta_phi_largest   = (delta_phi_largest - min(delta_phi_largest(:)))/(max(delta_phi_largest(:))-min(delta_phi_largest(:)));
        [dx,dy]=gradient(phi_single);
        mag = sqrt(dx.^2+dy.^2 + 0.001);
        curvature_largest = -divergence(dx./mag,dy./mag);
        if use_convex_hull
            kappa_weighted_edge_largest = (~CH_largest).*curvature_largest.*delta_phi_largest.*(curvature_largest>0);
        else
            kappa_weighted_edge_largest = curvature_largest.*delta_phi_largest.*(curvature_largest>0);
        end
        % Compute the SDF for the remaining components
        
        other_bwI = binaryImg - largest_bwI;
        CH_other = bwconvhull(other_bwI);
        phi_other =  double(bwdist(~other_bwI) - bwdist(other_bwI));
        delta_phi_other = Dirac(phi_other,epsilon);
%         delta_phi_other = imadjust(delta_phi_other,stretchlim(delta_phi_other),[0 ;1]);
        
        [dx,dy] = gradient(phi_other);
        mag = sqrt(dx.^2+dy.^2 + 0.001);
        curvature_other = -divergence(dx./mag,dy./mag);
        nx      = -dx.*delta_phi_other;
        ny      = -dy.*delta_phi_other;     % negative sign because LSF is positive inside, hence concave
        
         
        % Compute the attraction field using VFC
        
        rad = max(round(15*epsilon),1);
        Ker = AM_VFK(2,rad,'gaussian',gamma);   % Gaussian kernel of radii = rad
        
        vfkX = Ker(:,:,1);
        vfkY = Ker(:,:,2);
        
        Fx_master = conv2(kappa_weighted_edge_largest,vfkX,'same').*(~largest_bwI);
        Fy_master = conv2(kappa_weighted_edge_largest,vfkY,'same').*(~largest_bwI);
        
        % Compute for Slaves
        
        %     Fx_slaves = -conv2(kappa_weighted_edge_other,vfkX,'same').*(~other_bwI);
        %     Fy_slaves = -conv2(kappa_weighted_edge_other,vfkY,'same').*(~other_bwI);
        
        Fx_slaves = nx;
        Fy_slaves = ny;
        % Compute the attraction force
        
        attraction_force = (Fx_master.*Fx_slaves + Fy_master.*Fy_slaves).*curvature_other;
        
        %     imshow(kappa_weighted_edge_other,[]);drawnow;
        if max(attraction_force(:) ~=0)
            disp('attraction');
            attraction_force = (attraction_force - min(attraction_force(:)))/(max(attraction_force(:))-min(attraction_force(:)));
        end
        
        F_total = F_total + attraction_force;
    end
    %     attraction_force = attraction_force.*(delta_phi_other);
    %     attraction_force = attraction_force.*curvature_other;
    
    F_total = F_total/num_masters;
end

end

