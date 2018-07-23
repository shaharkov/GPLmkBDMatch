function H = ComputeDistHist(G, idx, num_bins, radius)

if ~exist('num_bins','var')
    num_bins = 20;
end
if ~exist('radius','var')
    radius = sqrt(G.ComputeSurfaceArea());
end

A = G.Triangulation2AdjacencyWeighted();


D = zeros(length(idx));
for ii = 1:length(idx)
    currDists = graphshortestpath(A,idx(ii));
    D(:,ii) = currDists(idx);
end

H = hist(D, linspace(0,radius,num_bins))';