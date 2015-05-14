function c = convergence(p_mask,n_mask,thresh,c)
    diff = p_mask - n_mask;
    n_diff = sum(abs(diff(:)));
    if n_diff < thresh
        c = c + 1;
    else
        c = 0;
    end
end

