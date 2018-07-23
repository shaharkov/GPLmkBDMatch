function [matchInds,regParam, isoParam] = matchSurfaceLmksBD(Gs,Lmks,putativeMatches,params)


%% parameters
visualize = params.visualize;
paramBDfilt = params.paramBDfilt;


%% setup


%% parameterize by minimizing isometric distortion:
% setup optimization problem
num_iter = 500;
TolX = 1e-10;
TolFun = 1e-6;
for j = 1:2
    % initialize with tutte
    [Gs{j}.Aux.V0, Gs{j}.Aux.inds_bF] = computeTutte(Gs{j}.F',Gs{j}.V',false); % global scale to be as ismetric as possible
    % constrain the centroid
    Gs{j}.Aux.eq_lhs_centroid = kron(eye(2),sparse(1,Gs{j}.Aux.inds_bF,1,1,size(Gs{j}.V,2)));
    Gs{j}.Aux.eq_rhs_centroid = [0;0];
    % setup optimization problem
    optimProblemParam{j} = OptimProblemIsoDist(Gs{j}.V', Gs{j}.F', Gs{j}.Aux.eq_lhs_centroid, Gs{j}.Aux.eq_rhs_centroid, Gs{j}.Aux.V0);
    % setup solver
    solverParam{j} = OptimSolverAcclQuadProx('QP', optimProblemParam{j}, false, true, true);
    % solve
    logParam{j} = solverParam{j}.solveTol(TolX, TolFun ,num_iter);
    % return isometric parameterization
    isoParam{j} = Mesh('VF',[solverParam{j}.X';zeros(1,Gs{j}.nV)],Gs{j}.F);
end

%% visualize putative correspondences
% if visualize
%     for j=1:2
%         putativeInds{j} = Lmks{j}(putativeMatches{j});
%     end
%     visualizeCorrepsondences2D(Gs, isoParam, putativeInds, false, 'Putative correpsondences')
% end


%% apply l0-feature matching with bounded distortion maps
if ~params.forceIdentity
% prepare
leftInd = colStack( Lmks{1}( putativeMatches{1} ) )';
leftBoundaryInd = Gs{1}.FindOrientedBoundaries;
rightInd = colStack( Lmks{2}( putativeMatches{2} ) )';
rightCoords = isoParam{2}.V(1:2, rightInd);
% retriangulate landmark points + boundary
leftV_ = isoParam{1}.V(1:2,[Lmks{1}, leftBoundaryInd]);
leftF_ = delaunay(leftV_')';
leftInd_ = putativeMatches{1}(:)';
leftBoundaryInd_ = length(Lmks{1}) + (1:length(leftBoundaryInd));
% BD filtering
matches = BDMatchMODIFIED(...
    leftF_, leftV_, leftInd_, leftBoundaryInd_, rightCoords,...
    paramBDfilt.K ,paramBDfilt.pnorm, paramBDfilt.initdelta, paramBDfilt.boundaryConstraintType, visualize);
matchInds{1} = leftInd(matches);
matchInds{2} = rightInd(matches);
else
    % override matching
    for j = 1:2
        matchInds{j} = Lmks{j};
    end
end


%% interpolate a map -- deform based on landmarks (minimize isometric distortion)
% global alignment -- apply global rigid tranformation to {1} based on correspondences
[T.R, T.t] = findBestRigidMotion(isoParam{1}.V(1:2,matchInds{1})',isoParam{2}.V(1:2,matchInds{2})',1);
% compute linear system (corresponding to landmark matches)
V = (isoParam{1}.V(1:2,:)'+T.t)*T.R';
if det(T.R)>0
    F =  isoParam{1}.F';
else
    F =  isoParam{1}.F([1 3 2],:)';
end
[eq_lhs,eq_rhs] = indCoordsToLinearSystem(isoParam{1}.V(1:2,:)',...
    matchInds{1}, ...
    isoParam{2}.V(1:2,matchInds{2})');
% setup optimization problem
num_iter = 500;
TolX = 1e-10;
TolFun = 1e-6;
optimProblem = OptimProblemIsoDist(V, F, eq_lhs, eq_rhs, [], 5);
% setup solver
solverInterp = OptimSolverAcclQuadProx('QP', optimProblem, false, true, true);
% solve
logInterp = solverInterp.solveTol(TolX, TolFun ,num_iter);
% store registered parameterizations
regParam{1} = Mesh('VF',[solverInterp.X';zeros(1,isoParam{1}.nV)],isoParam{1}.F);
regParam{2} = Mesh('VF',isoParam{2}.V,isoParam{2}.F);


%% plot parameterizations
if visualize
    % visualize
    %visualizeSolversForExamples;
    figure;
    patch('vertices', regParam{2}.V(1:2,:)', 'faces', regParam{2}.F', 'facecolor', [1,0,0], 'facealpha', 0.1, 'edgealpha', 0.1)
    hold on;
    patch('vertices', regParam{1}.V(1:2,:)', 'faces', regParam{1}.F', 'facecolor', [0,0,1], 'facealpha', 0.1, 'edgealpha', 0.1)
    scatter(regParam{1}.V(1,Lmks{1}), regParam{1}.V(2,Lmks{1}), 100, 'b'); %, 'filled')
    scatter(regParam{2}.V(1,Lmks{2}), regParam{2}.V(2,Lmks{2}), 'r'); %, 'filled')
    scatter(regParam{1}.V(1,matchInds{1}), regParam{1}.V(2,matchInds{1}), 90, 'b', 'filled')
    scatter(regParam{2}.V(1,matchInds{2}), regParam{2}.V(2,matchInds{2}), 30, 'r', 'filled')
    axis equal;
end