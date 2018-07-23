function [distBD, Lmks, putativeInds, matchIndsBD, regParamBD, isoParamBD] = wrapperLmkMatch(Gs,params)

% compute GP landmarks
for j = 1:2
    switch params.typeLmk
        case 'GP'
            Lmks{j} = Gs{j}.GetGPLmk(params.GPLmk.numLmk);
        case 'GP_nW'
            Lmks{j} = Gs{j}.GetGPLmk_NoCurvature(params.GPLmk.numLmk);
        case 'GP_Euc'
            Lmks{j} = Gs{j}.GetGPLmk_Euclidean(params.GPLmk.numLmk);
        case 'GT'
            Lmks{j} = Gs{j}.Aux.GTLmks(1:params.GTLmk.numLmk);
        otherwise
            error('invalid typeLmk')
    end
    Lmks{j} = reshape(Lmks{j},1,[]);
end
% compute putative correspondences
[putativeMatches, putativeInds] = computePutativeMatches(Gs,Lmks,params.computePutativeMatches);
% compute BD filtered landmark matches
[matchIndsBD,regParamBD,isoParamBD] = matchSurfaceLmksBD(Gs,Lmks,putativeMatches,params.matchSurfaceLmksBD);
% compute the Induced CP Distance
distBD = MapToDist(Gs{1}.V,Gs{2}.V,...
    knnsearch(regParamBD{2}.V(1:2,:)',regParamBD{1}.V(1:2,:)'),...
    Gs{1}.Aux.VertArea);
fprintf('Lmk-BD: Induced CP distance %f\n', distBD);