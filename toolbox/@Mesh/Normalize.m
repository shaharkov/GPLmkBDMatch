function Normalize(G)

% normalize
[G.Aux.Area,G.Aux.Center] = G.Centralize('ScaleArea');
% compute vertex areas
[~,TriArea] = G.ComputeSurfaceArea();
G.Aux.VertArea = (1/3)*G.F2V'*TriArea(:);

end