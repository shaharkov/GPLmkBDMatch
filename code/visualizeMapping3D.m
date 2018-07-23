function h = visualizeMapping3D(Gs, paramMeshes, refParam, titleText, params)

%% parameters
res = 500;
patch_option.SpecularStrength = .3;

%% setup texture
I = imread('texture/uv_grid2.jpg');
% I = fliplr(I);
I = imresize(I,[res res]);

% compute textures wrt reference parameterization
if isempty(refParam)
    Texture{1} = paramMeshes{1}.V(1:2,:)';
    Texture{2} = paramMeshes{2}.V(1:2,:)';
else
    Texture{1} = applyPiecewiseLinearMap(paramMeshes{2}.V(1:2,:)',refParam.V(1:2,:)',refParam.F',paramMeshes{1}.V(1:2,:)');
    Texture{2} = refParam.V(1:2,:)';
end

% scale
scale = max(sqrt(sum(Texture{2}.^2,2)));


%% draw
figure;
for j=1:2
    subplot(1,2,j);
    hS = approxTextureMap(Texture{j},Gs{j}.V',I,scale);
    hS.SpecularStrength = patch_option.SpecularStrength;
    axis off;
    axis equal;
    axis tight;
    axis xy;
    cameratoolbar;
    lighting phong
    camlight('headlight');
    camlight(180,0);
    h(j) = gca;
end
suptitle(titleText);
