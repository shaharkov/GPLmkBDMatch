%% init
clear
close all;
initialize;
yalmip('clear');
rng(1);


%% data
dataFolder = 'sample_data/';
params.GTLmk.numLmk = 16;
dataFileNames = {'a16_sas_aligned', 'b02_sas_aligned'};


%% Landmark Matching Method (see Table 1 in the paper)
params.typeLmkMatching = 'GP_BD'; % {'GP_BD','GP_Euc_BD','GP_nW_BD','GT_BD','GT2_BD','GT3'};


%% parameters
visualize = true;

params.preprocess.force = false;
params.preprocess.suffix = '_PRE';
params.preprocess.dataFileType = 'off';
params.preprocess.CPM.NumDensityPnts = 100;

params.CPM.FeatureType = 'ConfMax';
params.CPM.AngleIncrement = 0.05;
params.CPM.NumFeatureMatch = 4;
params.CPM.GaussMinMatch = 'off';

params.GPLmk.numLmk = 40;
params.computePutativeMatches.forceIdentity = false; % enabling makes sense only for GT landmarks
params.computePutativeMatches.WKSNN = 2;
params.matchSurfaceLmksBD.visualize = visualize;
params.matchSurfaceLmksBD.forceIdentity = false; % enabling makes sense only for GT landmarks
params.matchSurfaceLmksBD.paramBDfilt.K = 1.5;
params.matchSurfaceLmksBD.paramBDfilt.pnorm = 0.1;
params.matchSurfaceLmksBD.paramBDfilt.initdelta = 1;
params.matchSurfaceLmksBD.paramBDfilt.boundaryConstraintType = 'similarity+translation'; %{'none','fixed','linear','similarity','affine','similarity+translation'}

params.visualizeCorrepsondences2D.offsetFactors = [-1.1 0];
params.visualizeCorrepsondences2D.rotationAngle = 0*(pi/180);

% update parameters according to selected method
params = updateParamsPerMethod(params, params.typeLmkMatching);


%% load meshes and preprocess
for j = 1:2
     Gs{j} = loadAndPreprocess(dataFolder,dataFileNames{j},params.preprocess);
end


%% CPM
[distCPM, regParamCP, CPM12] = wrapperCPM(Gs,params);


%% Lmk Matching
[distBD, Lmks, putativeInds, matchIndsBD, regParamBD, isoParamBD] = wrapperLmkMatch(Gs,params);


%% ODLP
[distODLP] = wrapperODLP(Gs,params);



%% visualize results
% 2d parameterizations
visualizeCorrepsondences2D(Gs, isoParamBD, putativeInds, false, 'Putative correspondences (iso parameterization)', params.visualizeCorrepsondences2D);
visualizeCorrepsondences2D(Gs, regParamBD, putativeInds, false, 'Putative correspondences (registered parameterization)', params.visualizeCorrepsondences2D);
visualizeCorrepsondences2D(Gs, regParamBD, matchIndsBD, true, 'BD-filtered correpsondences', params.visualizeCorrepsondences2D);
% 3d and texture map
h{1} = visualizeCorrespondences3D(Gs, Lmks, false, 'Landmarks');
h{2} = visualizeCorrespondences3D(Gs, matchIndsBD, true, 'BD-filtered correpsondences');
h{3} = visualizeMapping3D(Gs, regParamBD, [], sprintf('Lmk-BD: Induced CP distance = %f', distBD));
h{4} = visualizeMapping3D(Gs, regParamCP, regParamBD{2}, sprintf('CPM: Induced CP distance = %f', distCPM));
% link
linkprop([h{:}], {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});