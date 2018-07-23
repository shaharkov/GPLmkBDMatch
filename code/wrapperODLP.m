function [distODLP] = wrapperODLP(Gs,params)

[~,~,distODLP] = findBestRigidMotion(Gs{1}.V(:,Gs{1}.Aux.GTLmks(1:params.GTLmk.numLmk))',...
    Gs{2}.V(:,Gs{2}.Aux.GTLmks(1:params.GTLmk.numLmk))', (1/params.GTLmk.numLmk)*ones(1,params.GTLmk.numLmk));
fprintf('ODLP distance %f\n', distODLP);
