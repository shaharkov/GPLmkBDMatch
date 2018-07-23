function [distCPM, regParamCP, CPM12] = wrapperCPM(Gs,params)

% compute CP maps
CPM12 = Gs{1}.ComputeContinuousProcrustes(Gs{2},params.CPM);
% extract texture maps
for j = 1:2
    regParamCP{j} = Mesh('VF',[CPM12.(['TextureCoords' num2str(j)]); zeros(1,Gs{j}.nV)], Gs{j}.F);
end
% compute the Induced CP Distance
distCPM = MapToDist(Gs{1}.V,Gs{2}.V,...
    knnsearch(regParamCP{2}.V(1:2,:)',regParamCP{1}.V(1:2,:)'),...
    Gs{1}.Aux.VertArea);
fprintf('CPM: Induced CP distance %f\n', distCPM);
fprintf('CPM distance %f\n', CPM12.cPdist);