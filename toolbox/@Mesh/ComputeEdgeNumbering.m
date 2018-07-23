function [M, edge_to_vertex_map, edge_num, M_orientation] = ComputeEdgeNumbering(G)
%use Gabriel Peyr's code
F = G.F;
Nv = size(G.V,2);

i = [F(1,:) F(2,:) F(3,:) F(2,:) F(3,:) F(1,:)];
j = [F(2,:) F(3,:) F(1,:) F(1,:) F(2,:) F(3,:)];
I = find(i<j);
i = i(I); j = j(I);
[~,I] = unique(i + 1234567*j);
i = i(I); j = j(I);
ne = length(i); % number of edges
s = (1:ne);

A = sparse([i;j],[j;i],[s;s],Nv,Nv);

edge_num = ne;
M=A;

edge_to_vertex_map = [i' j'];

%on the boundary make sure the edges are positively oriented
i = [F(1,:) F(2,:) F(3,:) ];
j = [F(2,:) F(3,:) F(1,:) ];
M_orientation = sparse(i,j,ones(1,length(i)),Nv,Nv);
M_orientation = M_orientation - M_orientation';
[row, col] = find(M_orientation>10^-8);
edge_to_vertex_map(A(M_orientation>10^-8),:) = [row col];

