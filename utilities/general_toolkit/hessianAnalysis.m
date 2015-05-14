function [val] = hessianAnalysis(H,smin,smax)
%HESSIANANALYSIS Compute the vesselness at the different scales for a voxel

    val = 0;
    [V,D] = eig(H);
    alpha = 0.5;
    beta = 0.5;
    c = 0.5;

    D = diag(D);
    [~,sortindex] = sort(abs(D));
    
    if D(sortindex(2))< 0 && D(sortindex(3)) < 0
        l1 = D(sortindex(1));
        l2 = D(sortindex(2));
        l3 = D(sortindex(3));
        
        R_A = abs(l2)/abs(l3);
        R_B = abs(l1)/sqrt(abs(l2)*abs(l3));
        S = sqrt(l1^2+l2^2+l3^2);
        
        val = (1 - exp(-(R_A^2)/(2*alpha^2)))*(exp(-(R_B^2)/(2*beta^2)))*(1-exp(-(S^2)/(2*c^2)));
        if isnan(val)
            fprintf('\n NAN');
            disp([R_A R_B S]);
        end
                
    end
    
end

