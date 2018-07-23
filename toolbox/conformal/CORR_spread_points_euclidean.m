function [sprd] = CORR_spread_points_euclidean(X,seed_inds,n)
%calc n farthest points with seed_inds as seed set and Euclidean distance
N = size(X,1);

dist = inf*ones(N,1);
if(~isempty(seed_inds) )

    for k=1:length(seed_inds)
        tdist = sum((X - ones(N,1)*X(seed_inds(k),:)).^2,2);
        dist = min([dist tdist],[],2);
    end
    
end

sprd = zeros(n,1);
for k=1:n
    progressbar(k,n,40);
    [~, tind] = max(dist);
    sprd(k) = tind(1);
    dist_from_k = sum((X - ones(N,1)*X(sprd(k),:)).^2,2);
    dist = min([dist dist_from_k],[],2);
end
    
    
end

