function [InterpInds1,InterpInds2,preInterpInds1] = FindEuclideanMutuallyNearestNeighbors(GM,GN,map,Type)
%FINDEUCLIDEANMUTUALLYNEARESTNEIGHBORS Summary of this function goes here
%   Detailed explanation goes here

switch Type
    case 'ADMax'
        GM_MaxInds = GM.Aux.ADMaxInds;
        GN_MaxInds = GN.Aux.ADMaxInds;
    case 'ConfMax'
        GM_MaxInds = GM.Aux.ConfMaxInds;
        GN_MaxInds = GN.Aux.ConfMaxInds;
    case 'GaussMax'
        GM_MaxInds = GM.Aux.GaussMaxInds;
        GN_MaxInds = GN.Aux.GaussMaxInds;
    case 'GaussMin'
        GM_MaxInds = GM.Aux.GaussMinInds;
        GN_MaxInds = GN.Aux.GaussMinInds;
end
pfGM_MaxInds = map(GM_MaxInds);

if ~isempty(GM_MaxInds)&&~isempty(GN_MaxInds)
    [~,R,~] = MapToDist(GM.V,GN.V,map,GM.Aux.VertArea);
    EucDistMatrix = pdist2((R*GM.V(:,GM_MaxInds))',GN.V(:,GN_MaxInds)')';
    
    [~, tind1] = min(EucDistMatrix,[],2);
    [~, tind2] = min(EucDistMatrix,[],1);
    tind2 = tind2';
    
    InterpMaxInds2 = find(tind2(tind1)==(1:length(tind1))');
    InterpMaxInds1 = tind1(InterpMaxInds2);
    
    InterpInds1 = pfGM_MaxInds(InterpMaxInds1);
    InterpInds2 = GN_MaxInds(InterpMaxInds2);
    preInterpInds1 = GM_MaxInds(InterpMaxInds1);
end

end

