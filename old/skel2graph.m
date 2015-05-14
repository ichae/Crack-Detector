function [skel_graph] = skel2graph(skel)
%SKEL2GRAPH Create a graph from a 3D skeleton
%           skel                    -- the binary skeleton
%           skel_graph.mat          -- non-sparse graph matrix created from the skeleton   
%           skel_graph.node_id      -- index of the nodes
%           skel_graph.node_coord   -- linear indices of each node coordinate

[n_r,n_c] = size(skel);
N = sum(skel(:));
G = zeros(N);
node_id = 1 : 1 : N;        % node id
node_coord = find(skel);    % coordinates of each node

for ii = 1 : N
        n1 = node_coord(ii);
        [r1,c1] = ind2sub(size(skel),n1);
        n1_S = sub2ind(size(skel),min(r1+1,n_r),c1);
        n1_N = sub2ind(size(skel),max(r1-1,1),c1);
        n1_E = sub2ind(size(skel),r1,min(c1+1,n_c));
        n1_W = sub2ind(size(skel),r1,max(c1-1,1));
        n1_NE = sub2ind(size(skel),max(r1-1,1),min(c1+1,n_c));
        n1_SE = sub2ind(size(skel),min(r1+1,n_r),min(c1+1,n_c));
        n1_NW = sub2ind(size(skel),max(r1-1,1),max(c1-1,1));
        n1_SW = sub2ind(size(skel),min(r1+1,n_r),max(c1-1,1));
    for jj = ii + 1 : N
        n2 = node_coord(jj);
        % Check for neighbours
        if (n2 == n1_S || n2 == n1_N || n2 == n1_E || n2 == n1_W || n2 == n1_NE || n2 == n1_NW || n2 == n1_SE || n2 == n1_SW)
            G(ii,jj) = 1;
        end
    end
end

G = G+G';

skel_graph.mat = G;
skel_graph.node_id = node_id;
skel_graph.node_coord = node_coord;


end

