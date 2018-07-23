function [VorArea] = ComputeVoronoiArea(G,Indices)
%COMPUTEVORONOIAREA Summary of this function goes here
%   Detailed explanation goes here

[~,~,Q] = G.PerformFastMarching(Indices);
VorArea = zeros(size(Indices));
for i=1:length(VorArea)
    VorArea(i) = sum(G.Aux.VertArea(Q == Indices(i)));
end

end

