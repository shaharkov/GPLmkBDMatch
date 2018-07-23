function [InterpInds1,InterpInds2,preInterpInds1] = FindMutuallyNearestNeighbors(GM,GN,map,Type)
%FINDMUTUALLYNEARESTNEIGHBORS Summary of this function goes here
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
    [~,~,Q] = GN.PerformFastMarching(pfGM_MaxInds);
    GN2pfGM = Q(GN_MaxInds);
    tind1 = zeros(size(GN_MaxInds));
    for j=1:length(tind1)
        tind1(j) = find(pfGM_MaxInds==GN2pfGM(j),1);
    end
    
    [~,~,Q] = GN.PerformFastMarching(GN_MaxInds);
    pfGM2GN = Q(pfGM_MaxInds);
    tind2 = zeros(size(pfGM_MaxInds));
    for j=1:length(tind2)
        if ~isempty(find(GN_MaxInds==pfGM2GN(j),1))
            tind2(j) = find(GN_MaxInds==pfGM2GN(j),1);
        else
            tind2(j) = NaN;
        end
    end
    
    InterpMaxInds2 = find(tind2(tind1)==(1:length(tind1))');
    InterpMaxInds1 = tind1(InterpMaxInds2);
    
    InterpInds1 = pfGM_MaxInds(InterpMaxInds1);
    InterpInds2 = GN_MaxInds(InterpMaxInds2);
    preInterpInds1 = GM_MaxInds(InterpMaxInds1);
end

end

