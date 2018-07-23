function E2F = ComputeE2F(obj)
% COMPUTEE2F computes E2F matrix: Ne by Nf matrix such that every entry
%             (e,f) = 1 iff  edge e is on face f, 0 otherwise.
%   input: TriangleMesh object
%   output: E2F as described above
%
% Created by Nave Zelzer on may 11 2014.
Nf = obj.nF;
Ne = obj.nE;
EFA = obj.ComputeEFA();

% now we get the edge indices for each face: the way we do it is not
% trivial. unique returns a vector of indices to the result unique vector
% thus for each edge (row) in unique(E(:,1:2)) this vector contains index
% this row
[~, ~, I] = unique(EFA(:,1:2), 'rows');
J = EFA(:,3);% third column contains sorted face indices
S = true(1,length(I));
E2F = sparse(I,J,S,Ne,Nf);