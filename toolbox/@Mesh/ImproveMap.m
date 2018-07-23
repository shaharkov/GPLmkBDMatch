function [rslt] = ImproveMap(GM,GN,DistMatrix,MapMatrix,TaxaCode,options)
%IMPROVEMAP
%   rslt.Gname1:            name of the first mesh
%   rslt.Gname2:            name of the second mesh
%   rslt.ImprDist:            continuous Procrustes distance
%   rslt.ImprMap:             optimal map generating cP distance
%   rslt.invImprMap:          inverse of rslt.ImprMap
%   rslt.TextureCoords1:    texture coordinates for the first mesh
%                           (deformed)
%   rslt.TextureCoords2:    textrue coordinates for the second mesh
%                           (not deformed)
%   rslt.ref:               =0 if Improved map is orientation-preserving
%                           =1 if Improved map is orientation-reversing
%
%   Tingran Gao, trgao10@math.duke.edu
%   last modified: 17 Aug 2014
%

if nargin<6
    options = [];
end
ProgressBar = getoptions(options,'ProgressBar','on');
ImprType = getoptions(options,'ImprType','MST');
SmoothMap = getoptions(options,'SmoothMap',0);
FeatureFix = getoptions(options,'FeatureFix','on');
if (SmoothMap==1) && (~isfield(options,'TextureCoords1Path')||~isfield(options,'TextureCoords2Path'))
    error('Need TextureCoords1Path and TextureCoords2Path provided in options fields!');
end
ChunkSize = getoptions(options,'ChunkSize',55);

if ~isfield(GM.Aux,'name') && ~isfield(GN.Aux,'name')
    error('Either Mesh missing .Aux.name');
end

rslt.Gname1 = GM.Aux.name;
rslt.Gname2 = GN.Aux.name;

TAXAind = cellfun(@(name) find(strcmpi(TaxaCode,name)),{GM.Aux.name,GN.Aux.name});
GroupSize = length(TaxaCode);

if ~strcmpi(ImprType,'Viterbi') && ~strcmpi(ImprType,'ComposedLAST') %%% MST or ComposedLAST
    ST = ConstructGraph(DistMatrix,ImprType,options);
    OptimalPath = FindGraphShortestPath(ST,TAXAind(1),TAXAind(2),TaxaCode,options);
    if SmoothMap==0
        %%% return vertex permutation map
        rslt.ImprMap = ComposeMapsAlongPath(OptimalPath,MapMatrix);
        [rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
        if det(R)>0
            rslt.ref = 0;
        else
            rslt.ref = 1;
        end
        rslt.TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
        rslt.TextureCoords2(:,isnan(compl(rslt.TextureCoords2))) = 1;
        if rslt.ref==1
            rslt.TextureCoords2(2,:) = -rslt.TextureCoords2(2,:);
        end
        rslt.TextureCoords1 = rslt.TextureCoords2(:,rslt.ImprMap);
    else
        %%% return a smooth map via texture coordinates interpolation
        [rslt.TextureCoords1,rslt.TextureCoords2] = ComposeTextureCoordsAlongPath(OptimalPath,options.TextureCoords1Path,options.TextureCoords2Path,GroupSize,ChunkSize,ProgressBar);
        rslt.ImprMap = knnsearch(rslt.TextureCoords2',rslt.TextureCoords1');
        [rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
        if det(R)>0
            rslt.ref = 0;
        else
            rslt.ref = 1;
        end
    end
elseif strcmpi(ImprType,'ComposedLAST') %%% cPLAST
    if ~exist(options.cPLASTPath,'file')
        ST = ConstructComposedLASTGraph(DistMatrix,MapMatrix,TaxaCode,options);
        save(options.cPLASTPath,'ST');
    else
        load(options.cPLASTPath);
    end
    OptimalPath = FindGraphShortestPath(ST,TAXAind(1),TAXAind(2),TaxaCode,options);
    if SmoothMap==0
        %%% return vertex permutation map
        rslt.ImprMap = ComposeMapsAlongPath(OptimalPath,MapMatrix);
        [rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
        if det(R)>0
            rslt.ref = 0;
        else
            rslt.ref = 1;
        end
        rslt.TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
        rslt.TextureCoords2(:,isnan(compl(rslt.TextureCoords2))) = 1;
        if rslt.ref==1
            rslt.TextureCoords2(2,:) = -rslt.TextureCoords2(2,:);
        end
        rslt.TextureCoords1 = rslt.TextureCoords2(:,rslt.ImprMap);
    else
        %%% return a smooth map via texture coordinates interpolation
        [rslt.TextureCoords1,rslt.TextureCoords2] = ComposeTextureCoordsAlongPath(OptimalPath,options.TextureCoords1Path,options.TextureCoords2Path,GroupSize,ChunkSize,ProgressBar);
        rslt.ImprMap = knnsearch(rslt.TextureCoords2',rslt.TextureCoords1');
        [rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
        if det(R)>0
            rslt.ref = 0;
        else
            rslt.ref = 1;
        end
    end
elseif strcmpi(ImprType,'Viterbi') %%% Viterbi
    ST = ConstructGraph(DistMatrix,'MST'); %%% MST provides a quick and dirty upper bound for Viterbi
    MSTPath = FindGraphShortestPath(ST,TAXAind(1),TAXAind(2),TaxaCode,options);
    MSTMap = ComposeMapsAlongPath(MSTPath,MapMatrix);
    [~,R] = MapToDist(GM.V,GN.V,MSTMap,GM.Aux.VertArea);
    options.R = R;
    options.T = length(MSTPath)-1;
    options.distMatrix = DistMatrix;
    options.mapMatrix = MapMatrix;
    options.SmoothMap = SmoothMap;
    [OptimalPath,rslt.ImprMap] = FindViterbiPath(GM,GN,DistMatrix,MapMatrix,TaxaCode,options);
    if SmoothMap==0
        %%% return vertex permutation map
        [rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
        if det(R)>0
            rslt.ref = 0;
        else
            rslt.ref = 1;
        end
        rslt.TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
        rslt.TextureCoords2(:,isnan(compl(rslt.TextureCoords2))) = 1;
        if rslt.ref==1
            rslt.TextureCoords2(2,:) = -rslt.TextureCoords2(2,:);
        end
        rslt.TextureCoords1 = rslt.TextureCoords2(:,rslt.ImprMap);
    else
        %%% return a smooth map via texture coordinates interpolation
        [rslt.TextureCoords1,rslt.TextureCoords2] = ComposeTextureCoordsAlongPath(OptimalPath,options.TextureCoords1Path,options.TextureCoords2Path,GroupSize,ChunkSize,ProgressBar);
        rslt.ImprMap = knnsearch(rslt.TextureCoords2',rslt.TextureCoords1');
        [rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
        if det(R)>0
            rslt.ref = 0;
        else
            rslt.ref = 1;
        end
    end
end

if strcmpi(FeatureFix,'on')
    disp('Performing Feature Fixing...');
    [rslt.TextureCoords1,rslt.ImprMap] = TPSDeformation(GM,GN,rslt.ImprMap,'ConfMax',rslt.TextureCoords1,rslt.TextureCoords2,options);
    disp('Done.');
end

rslt.invImprMap = knnsearch(rslt.TextureCoords1',rslt.TextureCoords2');

end

