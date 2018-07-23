function Ic = Adjacency2Incidence(G,A)

if nargin<2
    A = G.A;
end

%% compute list of edges
[i,j,s] = find(sparse(A));
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
Ic = sparse(is,js,s,n,nverts);
Ic = Ic';

% fix self-linking problem (0)
a = find(i==j);
if not(isempty(a))
    for t=a'
        Ic(i(t),t) = 1;
    end
end
