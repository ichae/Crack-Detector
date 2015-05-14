function [phi] = runTuFF2D(obj,param_Tuff)
%TUFF2DISOTROPIC Tubularity Flow Field in 2D, isotropic form
% $ \phi_t = (R(x)+\lambda \kappa +\nu F_{attr})|\nabla \phi|$

f = obj.Img;                                %-- actual smoothed image

phi0 = param_Tuff.phi0;                     %-- initialized sdf
R = param_Tuff.enhI;                        %-- the vesselness response map
smooth_trm = param_Tuff.smooth_trm ;        %-- smoothness coeff
attr_trm = param_Tuff.attr_trm ;            %-- attr. coeff
max_iter = param_Tuff.max_iter;             %-- iteration upperbound
error = param_Tuff.error;                   %-- convergence error
count_lim = param_Tuff.cnvg_cnt;            %-- convergence count
display = param_Tuff.interval;              %-- display interval
disp_Img = param_Tuff.disp_Img ;            %-- image to display 
magnify = param_Tuff.magnify ;              %-- display magnification
col = param_Tuff.col;                       %-- contour color



stop  = 0;
count = 0;
its   = 0;
phi   = phi0;


while (its < max_iter && stop == 0)
    
    idx = find(phi <= 1.5 & phi >= -1.5);   %-- get the curve's narrow band
    if ~isempty(idx)
        kappa = lsfCurvature(phi);
        [phix,phiy] = gradient(phi);
        grad_phi = sqrt(phix.^2+phiy.^2);
        
        dphi_dt1 = R.*grad_phi;
        dphi_dt2 = kappa.*grad_phi;
        
        dphi_dt1 = dphi_dt1/(eps+max(abs(dphi_dt1(:))));
        dphi_dt2 = dphi_dt2/(eps+max(abs(dphi_dt2(:))));
        
        if attr_trm == 0
            dphi_dt3 = 0;
        else
            bwI = phi>0;
            CC = bwconncomp(bwI);
            if CC.NumObjects > 1
                dphi_dt3 = TuFF.attrForce(bwI,param_Tuff.attr_range);
                dphi_dt3 = dphi_dt3/(eps+max(abs(dphi_dt3(:))));
            else
                dphi_dt3 = 0;
            end
        end
        
        F = dphi_dt1 + smooth_trm*dphi_dt2 + attr_trm*dphi_dt3;
        dt = 8/(eps+max(abs(F(:))));
        
        prev_mask = phi >=0 ;
        
        phi(idx) = phi(idx)+dt*F(idx);
        
        phi = TuFF.SussmanReinitLS(phi,0.5);
        phi = TuFF.NeumannBoundCond(phi);
        curr_mask = phi >=0 ;
        
        if display > 0
            if mod(its,display) == 0
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


function k = lsfCurvature(u)

f = fspecial('gaussian',6,2);
u = imfilter(u,f,'same');
[ux,uy] = gradient(u);
mag = eps+sqrt(ux.^2+uy.^2);
[uxx,uxy] = gradient(ux);
[uyx,uyy] = gradient(uy);

k = (uy.^2.*uxx+ux.^2.*uyy - ux.*uy.*uxy - ux.*uy.*uyx)./(mag).^3;
end

