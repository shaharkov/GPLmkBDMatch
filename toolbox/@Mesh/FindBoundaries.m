function [BV,BE] = FindBoundaries(obj)
% find the boundary vertices of G
% out: BV: indices of the boundary vertices
% out: BE: indices of the boundary edges
% Created by Nave Zelzer on may 22 2014.
EFA = obj.ComputeEFA();
% returns edge index for every edge in EFA, thus if we have duplicates the
% indices will have duplicates too.
[~, ~, I] = unique(EFA(:,1:2), 'rows');
% count how many of each edge
binc = histc(I,unique(I));
% if only one edge in the count we have a boundary.

F = obj.F;
% we pair up every triangle indices and create sparse matrix out of it
I = [F(1,:),F(2,:),F(3,:)];
J = [F(2,:),F(3,:),F(1,:)];
E = [I ; J];
% we sort every row such that (i,j) where i<j thus, E(i,1) < E(i,2)
E = sort(E)';
% sorting the rows lexicographically
E = sortrows(E);
% removing duplicates
obj.E = unique(E,'rows','stable')';

BE = obj.E(:,binc==1);
% BE = EFA(binc==1, 1:2)';
% get only the boundary vertices
BV = unique(BE(:))';
end