function EL = ComputeEdgeLengths(G)
V=G.V';
F=G.F';
L1 = sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2));
L2 = sqrt(sum((V(F(:,1),:)-V(F(:,3),:)).^2,2));
L3 = sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2));
EL = [L1, L2, L3];
end

