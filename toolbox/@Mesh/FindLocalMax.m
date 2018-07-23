function [ind,val] = FindLocalMax(G, f, K, exclude_boundary)
%Find local maxima of the scalar function f.

Nv = size(G.V,2);
BV = G.FindBoundaries;

AK=logical(G.A^K+speye(size(G.A))); % K-adjacency + self

max_f=max(AK*sparse(1:Nv,1:Nv,f),[],2);

if size(f,1)~=size(max_f,1)
    f = f';
end

if exclude_boundary==1
    ind = find(max_f==f);
    [i,j]=find(AK(BV,:)); % find vertices that are distant K from the boundary
    ind = setdiff(ind,j);
else
    ind = find(max_f==f);
end
val=f(ind);