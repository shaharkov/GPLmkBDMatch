function [LandmarkInds,Landmarks,GPNAS] = GetLandmarksFromPNAS(G,PathToLandmarks,PathToPNASMesh,NumLandmarks)
%GETLANDMARKS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    NumLandmarks = 16;
end

GPNAS = Mesh('off',PathToPNASMesh);
LandmarkFile = load(PathToLandmarks);
% tName = strtok(G.Aux.filename(end:-1:1), '/');
% tName = tName(end:-1:1);
% tName = strtok(tName, '_');
tName = G.Aux.name;
rawLandmarks = LandmarkFile.PP(strcmpi(LandmarkFile.names,tName),1:NumLandmarks,:);
Landmarks = zeros(size(rawLandmarks,2),3);
for k=1:size(rawLandmarks,2)
    Landmarks(k,:) = [rawLandmarks(1,k,1), rawLandmarks(1,k,2), rawLandmarks(1,k,3)];
end

[GPNAS.Aux.Area,GPNAS.Aux.Center] = GPNAS.Centralize('ScaleArea');
tree = KDTreeSearcher(GPNAS.V');
Landmarks = Landmarks-repmat(GPNAS.Aux.Center',NumLandmarks,1)*sqrt(1/GPNAS.Aux.Area);
LandmarkInds = tree.knnsearch(Landmarks);
Landmarks = GPNAS.V(:,LandmarkInds)';

end

