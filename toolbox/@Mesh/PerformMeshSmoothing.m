function f = PerformMeshSmoothing(G,f,options)

% perform_mesh_smoothing - smooth a function defined on a mesh by averaging
%
%   f = perform_mesh_smoothing(face,vertex,f,options);
%
%   Smooth a function f on a width of options.niter_averaging vertices.
%
%   Copyright (c) 2007 Gabriel Peyre

face = G.F;
vertex = G.V;

options.null = 0;
naver = getoptions(options, 'niter_averaging', 1);
% naver = 1;%getoptions(options, 'niter_averaging', 1);
type = 'combinatorial';%getoptions(options, 'averaging_type', 'combinatorial');
% type = 'cotangent'

if nargin<3
    f = [];
end
if isempty(f)
    f = vertex;
end
if size(f,1)<size(f,2)
    f = f';
end
[vertex,face] = G.CheckFaceVertex(vertex,face);

% if size(f,2)>1
%     for i=1:size(f,2)
%         f(:,i) = PerformMeshSmoothing(face,vertex,f(:,i),options);
%     end
%     return;
% end

n = max(face(:));

% compute normalized averaging matrix
if strcmp(type, 'combinatorial')
    %add diagonal
    W = G.Triangulation2Adjacency(face) + speye(n);
    D = spdiags(full(sum(W,2).^(-1)),0,n,n);
    W = D*W;
else
    %options.normalize=1;
    %W = compute_mesh_weight(vertex,face,type,options);
    W = G.Triangulation2AdjacencyWeighted;
    W = W + diag(sum(W));
    D = spdiags(full(sum(W,2).^(-1)),0,n,n);
    W = D*W;
end

% do averaging to smooth the field
for k=1:naver
    f = W*f;
end