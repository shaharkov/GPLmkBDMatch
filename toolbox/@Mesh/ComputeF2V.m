function F2V = ComputeF2V(G)

nf = size(G.F,2);
nv = size(G.V,2);
I = [G.F(1,:),G.F(2,:),G.F(3,:)];
J = [1:nf,1:nf,1:nf]';
S = ones(length(I),1);
F2V = sparse(J,I',S,nf,nv);
G.F2V = F2V;

end