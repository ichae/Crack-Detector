function [phi] = runFastEATuFF2D(self,param_Tuff)
%TUFF2DISOTROPIC Tubularity Flow Field with edge attraction in 2D, isotropic form
% Attraction term replaced by morphological closing

f = self.Img;                               %-- original smoothed image
phi0 = param_Tuff.phi0;                     %-- initialized sdf
R = param_Tuff.enhI;                        %-- the vesselness response map
smooth_trm = param_Tuff.smooth_trm ;        %-- smoothness coeff
edge_trm = param_Tuff.edge_trm ;            %-- edge. coeff
max_iter = param_Tuff.max_iter;             %-- iteration upperbound
error = param_Tuff.error;                   %-- convergence error
count_lim = param_Tuff.cnvg_cnt;            %-- convergence count
display = param_Tuff.interval;              %-- display interval
disp_Img = param_Tuff.disp_Img ;            %-- image to display 
magnify = param_Tuff.magnify ;              %-- display magnification
col = param_Tuff.col;                       %-- contour color

se = strel('disk',param_Tuff.attr_range);
% R = R/max(R(:));
%-- Edge Attraction Parameters
[fx,fy] = gradient(R);
edge_map = sqrt(fx.^2+fy.^2);
edge_map = edge_map/max(edge_map(:));       %-- edge map to compute edge vectors
rad = 30;
gamma = param_Tuff.edge_range;              %-- edge attraction range
Ker = TuFF.AM_VFK(2,rad,'gaussian',gamma);  % Gaussian kernel of radii = rad
vfkX = Ker(:,:,1); vfkY = Ker(:,:,2);
ex = conv2(edge_map,vfkX,'same');           %-- edge attraction vector
ey = conv2(edge_map,vfkY,'same');           %-- edge attraction vector

g = 1./(1+sqrt(fx.^2+fy.^2).^param_Tuff.p);
% g = g./max(g(:));
g_r = g.*(1-R);       % pulls the contour back from non vessel parts (high outside vessel)

stop  = 0;
count = 0;
its   = 0;
phi   = phi0;

restricted_pts = [];    
while (its < max_iter && stop == 0)
    
    idx = find(phi <= 1.5 & phi >= -1.5);   %-- get the curve's narrow band
    idx = setdiff(idx,restricted_pts);      %-- do not apply forces on these pixels
    if ~isempty(idx)
        kappa = lsfCurvature(phi);
        [phix,phiy] = gradient(phi);
        grad_phi = sqrt(phix.^2+phiy.^2);
        
        if param_Tuff.pull_back
            dphi_dt1 = (R-param_Tuff.pull_back*g_r).*grad_phi;            %-- TuFF force minus the g
            dphi_dt2 = g.*kappa.*grad_phi;      %-- Smoothness force
        else
            dphi_dt1 = (R).*grad_phi;           %-- TuFF force
            dphi_dt2 = kappa.*grad_phi;         %-- Smoothness force
        end
        
        dphi_dt1 = dphi_dt1/(eps+max(abs(dphi_dt1(:))));
        dphi_dt2 = dphi_dt2/(eps+max(abs(dphi_dt2(:))));
        
        dphi_dt4 = -(ex.*phix + ey.*phiy);  %-- Edge Attraction force
        dphi_dt4 = dphi_dt4/(eps+max(abs(dphi_dt4(:))));

        F = dphi_dt1 + smooth_trm*dphi_dt2 + edge_trm*dphi_dt4;
        dt = 5/(eps+max(abs(F(:))));
        
        prev_mask = (phi >=0) ;
        
        phi(idx) = phi(idx)+dt*F(idx);
        
        phi = TuFF.SussmanReinitLS(phi,0.5);
        phi = TuFF.NeumannBoundCond(phi);
        curr_mask = (phi >=0) ;
        
        if display > 0
            if mod(its,display) == 0
                [phi,pts] = morphAttraction(phi,se);
                restricted_pts = union(restricted_pts,pts);
                TuFF.showCurveAndPhi(phi,magnify,disp_Img,col);
                drawnow;
                fprintf('Iteration = %d\n',its);                
            end
        end
        
        count = TuFF.convergence(prev_mask,curr_mask,error,count);
        
        if count <= count_lim
            its = its + 1;
        else
            stop = 1;
            fprintf('Algorithm converged, iteration=%d \n',its);
        end
        
    else %-- Outside the narrow band
        break;
        
    end
    
end  %-- end of while loop

end


function [phi,restricted_pts] = morphAttraction(phi,se)

bwI = phi>=0;
clI = imclose(bwI,se);
d = (clI - bwI);
phi = TuFF.computeSDF(clI);
restricted_pts = find(d>0);

end


function k = lsfCurvature(u)

f = fspecial('gaussian',6,2);
u = imfilter(u,f,'same');
[ux,uy] = gradient(u);
mag = eps+sqrt(ux.^2+uy.^2);
[uxx,uxy] = gradient(ux);
[uyx,uyy] = gradient(uy);

k = (uy.^2.*uxx+ux.^2.*uyy - ux.*uy.*uxy - ux.*uy.*uyx)./(mag).^3;
end

