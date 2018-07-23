function oriented_BV = orient_boundary(BV,BE)
%ORIENT_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here

oriented_BV = zeros(size(BV));
BE(BE(:,1)==0,:) = [];
remaining = BV;
first_vertex = BV(1);
oriented_BV(1) = first_vertex;
remaining(remaining == first_vertex) = [];
for j=2:length(oriented_BV)
    cands = find(BE(:)==first_vertex);
    candRows = zeros(size(cands));
    candCols = zeros(size(cands));
    for k=1:length(cands)
        if (cands(k)>length(BV))
            candRows(k) = cands(k)-length(BV);
            candCols(k) = 2;
        else
            candRows(k) = cands(k);
            candCols(k) = 1;
        end
    end
    for k=1:length(cands)
        second_vertex = sum(BE(candRows(k),:))-first_vertex;
        if find(remaining == second_vertex)
            break;
        end
    end
    oriented_BV(j) = second_vertex;
    remaining(remaining == second_vertex) = [];
    first_vertex = second_vertex;
end


end

