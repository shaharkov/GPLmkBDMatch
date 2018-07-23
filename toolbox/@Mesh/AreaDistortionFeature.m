function AreaDistortions = AreaDistortionFeature(G,r)
%AREADISTORTIONFEATURE Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    BB = G.GetBoundingBox;
    diam = sqrt(sum(abs(BB(:,2)-BB(:,1)).^2));
    r = 0.05*diam;
end

if ~isfield(G.Aux,'VertArea')
    [~,TriArea] = G.ComputeSurfaceArea;
    G.Aux.VertArea = (TriArea'*G.F2V)'/3;
end

AreaDistortions = zeros(1,length(G.V));
for j=1:length(G.V)
    AreaDistortions(j) = sum(G.Aux.VertArea(find(sqrt(sum((G.V-repmat(G.V(:,j),1,length(G.V))).^2))<r)));
end

AreaDistortions = AreaDistortions/(pi*r^2);

end

