classdef TuFF < handle
    %TUFF
    % Code for enhancement via Frangi
    % Code for enhancement via LDE
    % Code for segmentation using TuFF (IEEE TIP,14)
    % Code for segmentation using EATuFF (IEEE Trans. ?)
    
    properties
        Img;   %-- original signal, not the vesselness
    end
    
    methods
        
        function obj = TuFF(I)
            obj.Img = I;
        end
        
        segI   = runTuFF2D(obj,param_Tuff)
        segI   = runEATuFF2D(obj,param_Tuff)
        segI   = runFastEATuFF2D(obj,param_Tuff)
        
    end % End methods
    
    
    methods(Static)
       
        phi0   = initRegionOtsu(Img,level,min_area)
        phi0   = initRegionManual(Img,type,magnify)

        enhImg = brightVesselEnhanceFrangi(I,scales, gamma)
        enhImg = darkVesselEnhanceFrangi(I,scales, gamma)
        enhImg = vesselEnhanceLDE(vess_param);

        F = attrForce(bwI,range)
        Ker = AM_VFK(dim,rad,type,gamma)
        
        function phi0 = computeSDF(bwI)
            % inside >= 0, outside < 0
            lsf = bwdist(bwI)-bwdist(1-bwI)+im2double(bwI)-.5;
            phi0 = -double(lsf); 
        end 
        
        
        function c = convergence(p_mask,n_mask,thresh,c)
            diff = p_mask - n_mask;
            n_diff = sum(abs(diff(:)));
            if n_diff < thresh
                c = c + 1;
            else
                c = 0;
            end
        end
        
        %-------------------------------------------------------------------------
        %               Check boundary condition
        function g = NeumannBoundCond(f)
            
            [nrow,ncol] = size(f);
            g = f;
            g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
            g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
            g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
        end
        
        function g = NeumannBoundCond3D(f)
            
            [nrow,ncol,depth] = size(f);
            g = f;
            g([1 nrow],[1 ncol],[1 depth]) = f([2 nrow-1],[2 ncol-1],[2 depth-1]);
            g(2:nrow-1,[1 ncol],[1 depth]) = f(2:nrow-1,[2 ncol-1],[2 depth-1]);
            g([1 nrow],2:ncol-1,[1 depth]) = f([2 nrow-1],2:ncol-1,[2 depth-1]);
            g([1 nrow],[1 ncol], 2:depth-1) = f([2 nrow-1],[2 ncol-1], 2:depth-1);
        end
        
        
        %-------------------------------------------------------------------------
        %          Reinitialize LSF by Sussman reinitialization method   
        
        function [D] = SussmanReinitLS(D,dt)
     
            %D  : level set function
            a = D - shiftR(D); % backward
            b = shiftL(D) - D; % forward
            c = D - shiftD(D); % backward
            d = shiftU(D) - D; % forward
            
            a_p = a;  a_n = a; % a+ and a-
            b_p = b;  b_n = b;
            c_p = c;  c_n = c;
            d_p = d;  d_n = d;
            
            a_p(a < 0) = 0;
            a_n(a > 0) = 0;
            b_p(b < 0) = 0;
            b_n(b > 0) = 0;
            c_p(c < 0) = 0;
            c_n(c > 0) = 0;
            d_p(d < 0) = 0;
            d_n(d > 0) = 0;
            
            dD = zeros(size(D));
            D_neg_ind = find(D < 0);
            D_pos_ind = find(D > 0);
            dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2)) - 1;
            dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)) - 1;
            
            D = D - dt .* TuFF.sussman_sign(D) .* dD;
        end
        
        
        function [D] = SussmanReinitLS3D(D,dt)
     
            %D  : level set function
            [nr,nc,nz] = size(D);
            a = D - cat(1, D(:,1,:), D(:,1:nc-1,:)); % backward
            b = cat(1, D(:,2:nc,:), D(:,nc,:)) - D;  % forward
            c = D - cat(2, D(1,:,:), D(1:nr-1,:,:)); % down
            d = cat(2, D(2:nr,:,:), D(nr,:,:)) - D; % down
            e = D - cat(3, D(:,:,1), D(:,:,nz-1));  % away- z dirn
            f = cat(3, D(:,:,2:nz), D(:,:,nz)) - D; % towards - z dirn
            
            a_p = a;  a_n = a; % a+ and a-
            b_p = b;  b_n = b;
            c_p = c;  c_n = c;
            d_p = d;  d_n = d;
            e_p = e;  e_n = n;
            f_p = f;  f_n = f;
            
            a_p(a < 0) = 0;
            a_n(a > 0) = 0;
            b_p(b < 0) = 0;
            b_n(b > 0) = 0;
            c_p(c < 0) = 0;
            c_n(c > 0) = 0;
            d_p(d < 0) = 0;
            d_n(d > 0) = 0;
            e_p(e < 0) = 0;
            e_n(e > 0) = 0;
            f_p(f < 0) = 0;
            f_n(f > 0) = 0;
            
            dD = zeros(size(D));
            D_neg_ind = find(D < 0);
            D_pos_ind = find(D > 0);
            dD(D_pos_ind) = sqrt(max(a_p(D_pos_ind).^2, b_n(D_pos_ind).^2) ...
                            + max(c_p(D_pos_ind).^2, d_n(D_pos_ind).^2) ...
                            + max(e_p(D_pos_ind).^2, f_n(D_pos_ind).^2)) - 1;
            dD(D_neg_ind) = sqrt(max(a_n(D_neg_ind).^2, b_p(D_neg_ind).^2) ...
                            + max(c_n(D_neg_ind).^2, d_p(D_neg_ind).^2)...
                            + max(e_n(D_neg_ind).^2, f_p(D_neg_ind).^2)) - 1;
            
            D = D - dt .* TuFF.sussman_sign(D) .* dD;
        end
        
        
        function S = sussman_sign(D)
            S = D ./ sqrt(D.^2 + 1);
        end
        
        %-------------------------------------------------------------------------
        %               Compute curvature of a function u(lsf)
        
        function k = compute_divergence(u)
            % Computes div(\nabla u|\nabla u|)
            [ux,uy] = gradient(u);
            normDu = sqrt(ux.^2+uy.^2+1e-10);	
            Nx = ux./normDu;
            Ny = uy./normDu;
            nxx = gradient(Nx);
            [~,nyy] = gradient(Ny);
            k = nxx+nyy;                        
        end
        
        
        %-------------------------------------------------------------------------
        %                   Display the evolving curve
        
        function showCurveAndPhi(phi,magnify,Img,cl)
            
            imshow(Img,[],'InitialMagnification',magnify);
            hold on;
            [c,h] = contour(phi,[0 0],cl,'Linewidth',3); hold off;
            test = isequal(size(c,2),0);
            while (test==false)
                s = c(2,1);
                if ( s == (size(c,2)-1) )
                    t = c;
                    hold on; plot(t(1,2:end)',t(2,2:end)',cl,'Linewidth',3);
                    test = true;
                else
                    t = c(:,2:s+1);
                    hold on; plot(t(1,1:end)',t(2,1:end)',cl,'Linewidth',3);
                    c = c(:,s+2:end);
                end
            end
        end
        
    end % End static methods
    
end % End class TuFF


%-------------------------------------------------------------------------
%                       Utility functions
%-------------------------------------------------------------------------

function shift = shiftD(M)
shift = shiftR(M')';
end
%-------------------------------------------------------------------------

function shift = shiftL(M)
shift = [ M(:,2:size(M,2)) M(:,size(M,2)) ];
end
%-------------------------------------------------------------------------

function shift = shiftR(M)
shift = [ M(:,1) M(:,1:size(M,2)-1) ];
end
%-------------------------------------------------------------------------

function shift = shiftU(M)
shift = shiftL(M')';
end
%--------



