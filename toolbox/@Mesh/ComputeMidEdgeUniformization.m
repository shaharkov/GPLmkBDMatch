function [uniV,uniF] = ComputeMidEdgeUniformization(G,options)
%COMPUTEMIDEDGEUNIFORMIZATION Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    options = [];
end

% SmoothCurvatureFields = getoptions(options,'SmoothCurvatureFields',10);
% DensityLocalWidth = getoptions(options,'DensityLocalWidth',5);
% ExcludeBoundary = getoptions(options,'ExcludeBoundary',1);

% [~,TriAreas] = G.ComputeSurfaceArea;
% G.Aux.VertArea = (TriAreas'*G.F2V)/3;
%%% compute mid-edge mesh
[mV,mF,M,E2Vmap] = G.ComputeMidEdge;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 1) decide where to cut the surface (for flattening/uniformization)
disp('Find a face to cut the surface for uniformization...');
[v_max_V] = CORR_spread_points_euclidean(G.V',[],200);
GeoField = pdist2(G.V(:,v_max_V)',G.V');

medianGeoField = mean(GeoField,2);
[~, minplc] = min(medianGeoField);
cut_vertex = v_max_V(minplc);
cut_face = find(G.F2V(:,cut_vertex),1);
disp('Found.');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 2) flatten the mesh conformally 
disp('Flattening the mid-edge mesh...')
unmV = CORR_map_mesh_to_plane_nonconforming(G.V',G.F',mF',cut_face,M,E2Vmap,G.nE,0);

unmF = mF';
unmF(cut_face,:) = []; %% it is the same face number as the original mesh

center_ind = cut_vertex;
tind = find(G.F2V(:,center_ind),1);%v_max_V1(kk);
center_ind = mF(1,tind);
disp('Flattened.');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 3) map domain to disk (add the infinity point back as sample point)
% transfer the indices of center point to the mid-edge
unmV = CORR_transform_to_disk_new_with_dijkstra(unmV,mF,E2Vmap,center_ind);
unmF = mF';
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 4) map the original mesh to the disk using the mid-edge structure
disp('Flattening the ORIGINAL mesh using the mid-edge flattening...')
uniV = CORR_flatten_mesh_by_COT_using_midedge(G.V',G.F',M,mV,unmF,unmV,cut_face);
uniF = G.F';
G.Aux.UniformizationV = [uniV,zeros(G.nV,1)]';
disp('Flattend.');
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 5) calculate conformal factor on G
Lambda = CORR_calculate_conformal_factors(G.F',G.V',uniV);
G.Aux.Conf = CORR_calculate_conformal_factors_verts(G.F',Lambda);
%%% fix nan in G.Aux.Conf
nind = find(isnan(G.Aux.Conf)); %all NaN vertices
if(~isempty(nind))
    PP = G.V'; %look for closest me-vert on original mesh which is not NaN later
    PP(nind,:)=[]; %remove the ones with NaN
    LL = G.Aux.Conf;
    LL(nind)=[];
    idxs = knnsearch(PP,G.V(:,nind)');
    G.Aux.Conf(nind) = LL(idxs);
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 6) spread density points on mesh
disp('Spreading density points on mesh...');
G.ExtractFeatures(options);
% Cgauss = G.ExtractFeatures(options);
% Conf = G.Aux.Conf;
% for j=1:SmoothCurvatureFields
%     WeightMatrix = repmat(TriAreas,1,G.nV).*G.F2V.*repmat(1./G.Aux.VertArea,G.nF,1);
%     
%     CgaussFace = mean(Cgauss(G.F));
%     Cgauss = CgaussFace*WeightMatrix;
%     
%     ConfFace = mean(Conf(G.F));
%     Conf = ConfFace*WeightMatrix;
% end
% [G.Aux.GaussMaxInds,~] = G.FindLocalMax(Cgauss',DensityLocalWidth,ExcludeBoundary);
% [G.Aux.GaussMinInds,~] = G.FindLocalMax(-Cgauss',DensityLocalWidth,ExcludeBoundary);
% [G.Aux.ConfMaxInds,~] = G.FindLocalMax(G.Aux.Conf',DensityLocalWidth,ExcludeBoundary);
minds = [G.Aux.GaussMaxInds;G.Aux.GaussMinInds;G.Aux.ConfMaxInds];
minds = unique(minds);

G.Aux.DensityPnts = G.GeodesicFarthestPointSampling(options.NumDensityPnts,minds);

disp('DONE!');

end

