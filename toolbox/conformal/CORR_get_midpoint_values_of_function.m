function [mu] = CORR_get_midpoint_values_of_function(F,u, M, e2v, nume)
%takes a function on a mesh and returns its values on the midedges (for the
%midedge mesh)

%[M e2v nume] = compute_edge_numbering(F);

mu = zeros(nume,1);

for k=1:nume
   vind=e2v(k,:);
    mu(k) = 0.5*(u(vind(1),:) + u(vind(2),:));
end 
    

