function visualizeCorrepsondences2D(Gs, paramMeshes, matchInd, colorVerts, titleText, params)

%% parameters
Ncam = [1 0 0];
patch_option = [];
patch_option.edgecolor = [0.1 0.1 0.1];
patch_option.edgealpha = 0.2;
patch_option.facecolor = 'flat';
if exist('params','var')
    offsetFactors = params.offsetFactors;
    rotationAngle = params.rotationAngle;
else
    offsetFactors = [-1.3 0];
    rotationAngle = 0;
end
R = blkdiag([cos(rotationAngle) sin(rotationAngle); -sin(rotationAngle) cos(rotationAngle)],1);
    

%% colors
if colorVerts
    c = repmat(jet(numel(matchInd{1})),[2 1]);
else
    c = 'b';
end

%% rotate meshes
for j = 1:2
    rotParamMeshes{j}.V = R*paramMeshes{j}.V;
end

%% offset meshes
visualizeOffset = [0 0 0;...
    [max(rotParamMeshes{2}.V(1,:))-min(rotParamMeshes{2}.V(1,:)), max(rotParamMeshes{2}.V(2,:))-min(rotParamMeshes{2}.V(2,:))].*offsetFactors, 0];
for j = 1:2
    offestParamMeshes{j} = Mesh('VF',rotParamMeshes{j}.V + visualizeOffset(j,:)',paramMeshes{j}.F);
end

%% plot
figure;
for j=1:2
    % normal map
    Nv = Gs{j}.ComputeNormal;
    patch_option.FaceVertexCData = 0.7 + 0.4*max((Ncam*Nv)',0)*[1 1 1];
    offestParamMeshes{j}.draw(patch_option);
end
hold on;
plot([offestParamMeshes{1}.V(1,matchInd{1}); offestParamMeshes{2}.V(1,matchInd{2})],...
    [offestParamMeshes{1}.V(2,matchInd{1}); offestParamMeshes{2}.V(2,matchInd{2})], 'linewidth',1, 'color',[0.5 0 0.4]);
scatter([offestParamMeshes{1}.V(1,matchInd{1}), offestParamMeshes{2}.V(1,matchInd{2})],...
    [offestParamMeshes{1}.V(2,matchInd{1}), offestParamMeshes{2}.V(2,matchInd{2})], 60, c, 'filled');
suptitle(titleText);