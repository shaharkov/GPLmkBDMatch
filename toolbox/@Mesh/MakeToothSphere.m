function sTooth = MakeToothSphere(G, numVertOnLid, numIntPts, shellRadius)
%MAKETOOTHSPHERE: close disk-type meshes, output the spherical version
%   numVertOnLid: number of points on the bottom
%   numIntPts:    number of (initial) points randomly sampled from the disk
%   shellRadius:  radius threshold for interior points in the disk (a real
%                 number between 0 and 1)
%
%   Tingran Gao (trgao10@math.duke.edu)
%   last modified: Sep 22, 2016
%

[oBV, oBE] = G.FindOrientedBoundaries();
arcLengths = cumsum([0,sqrt(sum(G.V(:, oBV(2:end)) - G.V(:, oBV(1:(end-1)))).^2)]);
totalLength = arcLengths(end) + norm(G.V(:,oBV(end))-G.V(:,oBV(1)));

%%%% preserve the orientation of the boundary when projected to the circle
pjV = G.V;
pjV(end, :) = 0;
Nv = cross(diff(pjV(:,oBE(:,1)), 1, 2), diff(pjV(:,oBE(:,2)), 1, 2));
if Nv(end) > 0
    [x,y] = pol2cart(arcLengths/totalLength*(2*pi), 1);
else
    [x,y] = pol2cart(-arcLengths/totalLength*(2*pi), 1);
end

intPts = rand([2,numIntPts]);
intPts = intPts * 2 - repmat(ones(2,1), 1, numIntPts);
intPts = intPts(:,sum(intPts.^2)<=shellRadius^2);

allPts = [[x;y],intPts];
DT = delaunayTriangulation(allPts');

discMesh = Mesh('VF', allPts, DT.ConnectivityList');

L = discMesh.ComputeCotanLaplacian();
b = zeros(size(allPts, 2), 3);
b(1:length(oBV), :) = G.V(:, oBV)';

L(1:length(oBV),:) = 0;
L(1:length(oBV), 1:length(oBV)) = eye(length(oBV));

reconCoords = L\b;

bMesh = Mesh('VF', reconCoords', discMesh.F);
VertSampInd = bMesh.GeodesicFarthestPointSampling(numVertOnLid, 1:length(oBV));
DT = delaunayTriangulation(discMesh.V(1:2, VertSampInd)');

%%%% adjust the orientation of the lid to stay consistent with respect to
%%%% the entire sphere-type mesh
LidFaceList = DT.ConnectivityList;
% checkOrienMesh = Mesh('VF', bMesh.V(:, VertSampInd'), DT.ConnectivityList);
% [~, NFs] = checkOrienMesh.ComputeNormal();
% NF = mean(NFs, 2);
% [~, origNFs] = G.ComputeNormal();
% origNF = mean(origNFs, 2);
% if NF(end)*origNF(end) > 0
%     LidFaceList = LidFaceList(:, [2 1 3]);
% end

% pjLidBV = bMesh.V(:, 1:length(oBV));
% pjLidBV(end,:) = 0;
% % Nv = cross(diff(pjLidBV(:, bMesh.F([1,2],1)), 1, 2), diff(pjLidBV(:, bMesh.F([2,3],1)), 1, 2));
% % Nv = cross([DT.Points(LidFaceList(1,2),:)-DT.Points(LidFaceList(1,1),:),0],...
% %            [DT.Points(LidFaceList(1,3),:)-DT.Points(LidFaceList(1,2),:),0]);
% Nv = cross(diff(pjLidBV(:, [10,20]), 1, 2), diff(pjLidBV(:, [20,30]), 1, 2));
% if Nv(end) > 0
%     LidFaceList = LidFaceList(:, [2 1 3]);
% end

% LidFaceList = LidFaceList(:, [1 3 2]);

% tSortLidFaceList = sort(LidFaceList, 2);
% orienFlagIdx = 0;
% for kk=1:size(tSortLidFaceList,1)
%     if (tSortLidFaceList(kk,1) <= length(oBV)) && (tSortLidFaceList(kk,2) <= length(oBV))
%         orienFlagIdx = kk;
%         break
%     end
% end
% orienPairInFlagRow = find(LidFaceList(orienFlagIdx, :) < length(oBV));
% orienPairInFlagRow = orienPairInFlagRow(1:2);
% if (orienPairInFlagRow(1) == 1) && (orienPairInFlagRow(2) == 3)
%     orienPairInFlagRow = [3,1];
% end
% orienFlagEdge = oBV(LidFaceList(orienFlagIdx, orienPairInFlagRow));
% if ~isempty(find(sum((oBE - repmat(orienFlagEdge',1,size(oBE,2))).^2)==0, 1))
%     LidFaceList = LidFaceList(:, [1 3 2]);
% end

% Nv = cross([DT.Points(LidFaceList(1,2),:)-DT.Points(LidFaceList(1,1),:),0],...
%            [DT.Points(LidFaceList(1,3),:)-DT.Points(LidFaceList(1,2),:),0]);
% if Nv(end) > 0
%     LidFaceList = LidFaceList(:, [2 1 3]);
% end

% transIdxTable = (1:numVertOnLid)+G.nV;
% transIdxTable(1:length(oBV)) = oBV;
% sTooth = Mesh('VF', [G.V, bMesh.V(:, VertSampInd)], [G.F, LidFaceList'+G.nV]);
% sTooth = Mesh('VF', [G.V, bMesh.V(:, VertSampInd)], [G.F, transIdxTable(LidFaceList')]);
transIdxTable = [oBV,G.nV+(1:numVertOnLid-length(oBV))];
sTooth = Mesh('VF', [G.V, bMesh.V(:, VertSampInd((length(oBV)+1):end))],...
                    [G.F, transIdxTable(LidFaceList')]);
% sTooth.F = sTooth.PerformFacesReorientation(sTooth.V, sTooth.F, struct('method','slow'));

%%%% pick an arbitrary boundary edge to test orientation compatiability
randomTestList = randi([1,size(oBE,2)], 1, 10);
compatibleYes = 0;
compatibleNo = 0;
for jj = 1:length(randomTestList)
    examEdge = oBE(:, randomTestList(jj));
    sortedEdgeList = sort(sTooth.E);
    examEdgeIdx = find(sum((sortedEdgeList-repmat(examEdge, 1, size(sortedEdgeList,2))).^2) == 0);
    examEdgeAdjFaceIdx = find(sTooth.E2F(examEdgeIdx, :));
    Face1 = sTooth.F(:, examEdgeAdjFaceIdx(1));
    Face2 = sTooth.F(:, examEdgeAdjFaceIdx(2));
    anchorIdx1 = find(Face1 == examEdge(1));
    anchorIdx2 = find(Face2 == examEdge(1));
    shiftAmount = anchorIdx2 - anchorIdx1;
    if shiftAmount < 0
        shiftAmount = shiftAmount + 3;
    end
    Face1 = circshift(Face1, [shiftAmount, 0]);
    if find(Face1 == examEdge(2)) == find(Face2 == examEdge(2))
        compatibleNo = compatibleNo + 1;
    else
        compatibleYes = compatibleYes + 1;
    end
%     examEdgeInFace1 = Face1(Face1 ~= setdiff(Face1, examEdge));
%     if (examEdge(1) == Face1(3)) && (examEdge(2) == Face1(1))
%         examEdgeInFace1 = examEdgeInFace1(end:-1:1);
%     end
%     examEdgeInFace2 = Face2(Face2 ~= setdiff(Face2, examEdge));
%     if (examEdge(1) == Face2(3)) && (examEdge(2) == Face2(1))
%         examEdgeInFace2 = examEdgeInFace2(end:-1:1);
%     end
%     if all(examEdgeInFace1 == examEdgeInFace2)
%         compatibleNo = compatibleNo + 1;
%     else
%         compatibleYes = compatibleYes + 1;
% %         sTooth.F(:, (G.nF+1):end) = sTooth.F([2 1 3], (G.nF+1):end);
%     end
end

if compatibleYes < compatibleNo
    sTooth.F(:, (G.nF+1):end) = sTooth.F([2 1 3], (G.nF+1):end);
end

end

