function [A,E] = ComputeAdjacencyMatrix(G)
if ~isempty(G.A)
    A=G.A;
end

V=G.V;
F=G.F;
Nv=size(V,2);
Nf=size(F,2);

I = [F(1,:),F(2,:),F(3,:)];
J = [F(2,:),F(3,:),F(1,:)];
S = [1:Nf,1:Nf,1:Nf];
E = sparse(I,J,S',Nv,Nv);
A = double(E > 0);
A = A + A';
A = double(A > 0);

G.A=A;
G.E=E;