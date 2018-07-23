function vring = ComputeVertexRing(G,face)

% compute_vertex_ring - compute the 1 ring of each vertex in a triangulation.
%
%   vring = compute_vertex_ring(face);
%
%   vring{i} is the set of vertices that are adjacent
%   to vertex i.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    face = G.F;
end
[~,face] = G.CheckFaceVertex([],face);

nverts = max(max(face));

A = G.Triangulation2Adjacency(face);
[i,j,s] = find(sparse(A));

% create empty cell array
vring{nverts} = [];

for m = 1:length(i)
    vring{i(m)}(end+1) = j(m);
end