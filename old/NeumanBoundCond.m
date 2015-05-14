function g = NeumanBoundCond(f)
    
    [nrow,ncol] = size(f);
    g = f;
%     g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
%     g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
%     g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

    g(1,1) = f(2,2); g(1,ncol) = f(2,ncol-1);
    g(nrow,1) = f(nrow-1,2); g(nrow,ncol) = f(nrow-1,ncol-1);
    
    g(1,2:ncol-1) = f(2,2:ncol-1);
    g(nrow,2:ncol-1) = f(nrow-1,2:ncol-1);
    
    g(2:nrow-1,1) = f(2:nrow-1,2);
    g(2:nrow-1,ncol) = f(2:nrow-1,ncol-1);
end
