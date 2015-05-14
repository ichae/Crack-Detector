function [n_t1,n_t2,n_reg] = normalizeTerms(t1,t2,t_reg,scale )
%NORMALIZETERMS normalize the terms

epsilon = 0.0001;

m1 = min(t1(:));
m2 = min(t2(:));
m3 = min(t_reg(:));

m = min([m1;m2;m3]);


M1 = max(t1(:));
M2 = max(t2(:));
M3 = max(t_reg(:));

M = max([M1;M2;M3]);


n_t1 = -scale + 2*((t1-m)/(M-m))*scale;
n_t2 = -scale + 2*((t2-m)/(M-m))*scale;
n_reg = -scale + 2*((t_reg-m)/(M-m))*scale;

end

