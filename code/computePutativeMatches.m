function [putativeMatches, putativeInds] = computePutativeMatches(Gs,Lmks,params)

if ~params.forceIdentity
    
    % compare WKS
    putativeMatches{2} = knnsearch(Gs{2}.Aux.WKS(Lmks{2},:), Gs{1}.Aux.WKS(Lmks{1},:), 'k', params.WKSNN);
    putativeMatches{1} = ndgrid(1:size(putativeMatches{2},1), 1:size(putativeMatches{2},2));
    for j=1:2
        putativeInds{j} = Lmks{j}(putativeMatches{j});
    end
else
    assert(isequal(size(Lmks{1}),size(Lmks{2})), 'Lmks must be of same size');
    for j=1:2
        putativeMatches{j} = (1:length(Lmks{1}))';
        putativeInds{j} = reshape(Lmks{j},[],1);
    end
end

