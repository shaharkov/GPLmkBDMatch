function h = visualizeCorrespondences3D(Gs, matchInd, colorVerts, titleText, params)

%% parameters
patch_option = '';
patch_option.edgecolor = 'none';
% patch_option.edgealpha = 0.1;
% patch_option.facealpha = 0.8;
patch_option.facecolor = [.9 .9 .8];
patch_option.SpecularStrength = .3;
sphereSize = 0.02;


%% colors
if colorVerts
    c = repmat(jet(numel(matchInd{1})),[2 1]);
else
    c = [0 0 1];
end


%% draw
figure;
for j=1:2
    subplot(1,2,j);
    Gs{j}.draw(patch_option);
    hold on
    drawSpheres(Gs{j}.V(:,matchInd{j})',sphereSize,c);
    lighting phong;
    camlight('headlight');
    camlight(180,0);
    h(j) = gca;
end
suptitle(titleText);
