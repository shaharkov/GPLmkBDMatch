function p_out = applyPiecewiseLinearMap(V1,V2,F,p_in)
% apply the mapping that takes the triangulation (V1,F) to (V2,F) on the
% point cloud p_in (2D)

Tr = triangulation(F,V1);
[ti,B] = Tr.pointLocation(p_in);

p_out = zeros(size(p_in));
for ii = 1:size(p_out,1)
    if ~isnan(ti(ii))
        p_out(ii,:) = B(ii,:)*V2(F(ti(ii),:),:);
    else
        p_out(ii,:) = -10;
    end
end