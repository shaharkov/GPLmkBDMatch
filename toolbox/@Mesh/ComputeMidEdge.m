function [mV,mF,M,e2v_map] = ComputeMidEdge(G)
%COMPUTEMIDEDGE Summary of this function goes here
%   Detailed explanation goes here

V = G.V;
F = G.F;
[M,e2v_map,~,~] = G.ComputeEdgeNumbering;

nume = size(G.V2E,2);
numf = size(G.F,2);

%the new vertices are the old edges
mV = zeros(3,nume);
for k = 1:nume
    vind = e2v_map(k,:);
    mV(:,k) = 0.5*(V(:,vind(1)) + V(:,vind(2)));
end

%the new faces are one per old face
mF = zeros(3,numf);
for k=1:numf
   f = F(:,k);
   
   ie1 = M(f(1),f(2));
   ie2 = M(f(2),f(3));
   ie3 = M(f(3),f(1));
    
   mF(:,k) = [ie1 ie2 ie3]';
end

end

