classdef DeShadow < handle
    %DESHADOW : Remove Shadow and denoise by variational method
    properties
        inImg;
        shadow;
    end
    
    methods
        
        % =========== CONSTRUCTOR ===============
        function obj = DeShadow(I)
            obj.inImg = I;
            obj.shadow = zeros(size(I));
        end
        
        % ========== Function Prototype =========
        
        function [s] = estimateShadowGaussian(obj,sigma)
            f = obj.inImg;
            sz = 6*sigma;
            h = fspecial('gaussian',sz,sigma);
            s = imfilter(f,h,'same','replicate');
        end
        s = estimateShadowLegendre(obj,n_poly)
        s = estimateShadowMinimax(obj,deshadow_param)
        s = reconstructImage(obj,deshadow_param);
        
        
    end
    
       % ========== Static Methods ===============
    methods(Static)
        
        function count = convergence(J,its,count,error)
            m = J(its);
            n = J(its-1);
            if abs(m-n) <= error
                count = count + 1;
            else
                count = 0;
            end
        end
        
        function g = NeumannBoundCond(f)
            [nrow,ncol] = size(f);
            g = f;
            g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
            g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
            g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
        end
        
    end
    
end

