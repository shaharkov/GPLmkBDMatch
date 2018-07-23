function EFA = ComputeEFA( obj )
% COMPUTEEFA compute an array that contains columns of the shape [i,j,f]
% where i,j is the indices of an edge and f is the face incident to that
% edge for all edges in the mesh.
%   input: TriangleMesh object
%   output: EFA as described above (EFA stands for edge face array)
%
% Created by Nave Zelzer on may 11 2014.
F = obj.F;
Nf = obj.nF;
% we pair up every triangle indices and create sparse matrix out of it
I = [F(1,:),F(2,:),F(3,:)];
J = [F(2,:),F(3,:),F(1,:)];
EFA = [I ; J];
% we sort every row such that (i,j) where i<j thus, E(i,1) < E(i,2)
EFA = sort(EFA);
% add face indices column (a third column)
EFA = [EFA ; 1:Nf,1:Nf,1:Nf]';
% sorting the rows lexicographically
EFA = sortrows(EFA);
end

