function V2E = ComputeV2E(G)
% Based on Gabriel Peyre's code
% imrpove - compute V2E directly
A=G.ComputeAdjacencyMatrix;
%% compute list of edges
[i,j] = find(sparse(A));
I = find(i>=j);
i = i(I);
j = j(I);
% number of edges
n = length(i);
% number of vertices
nverts = size(A,1);

%% build sparse matrix
s = [ones(n,1); -ones(n,1)];
is = [(1:n)'; (1:n)'];
js = [i(:); j(:)];
V2E = sparse(is,js,s,n,nverts);
V2E = V2E';

% fix self-linking problem (0)
a = find(i==j);
if not(isempty(a))
    for t=a'
        V2E(i(t),t) = 1;
    end
end
G.V2E=V2E;
end

