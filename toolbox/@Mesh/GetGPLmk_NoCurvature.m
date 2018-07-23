function [GPLmkIdx,ptuq] = GetGPLmk_NoCurvature(G,numLmk)
%GETGPLMK Summary of this function goes here
%   Detailed explanation goes here

% if nargin < 3
%     lambda = 0.5;
% end

% G.Centralize('ScaleArea');
% [~,TriArea] = G.ComputeSurfaceArea();
% G.Aux.VertArea = G.F2V'*TriArea;

% [Cgauss,Cmean] = G.ComputeCurvature();
% Lambda = G.Aux.VertArea.*(lambda*abs(Cgauss)/sum(abs(Cgauss))+(1-lambda)*abs(Cmean)/sum(abs(Cmean)));

%[~,curvature] = findPointNormals(G.V',10);
%Lambda = G.Aux.VertArea.*curvature/sum(curvature);
Lambda = G.Aux.VertArea;

if size(G.E,1) > 2
    [I,J] = find(tril(G.E));
    G.E = ([I,J])';
end
EdgeIdxI = G.E(1,:);
EdgeIdxJ = G.E(2,:);
bandwidth = mean(sqrt(sum((G.V(:,EdgeIdxI)-G.V(:,EdgeIdxJ)).^2)))/5;

BNN = min(500,G.nV);
% atria = nn_prepare(G.V');
% [idx, dist] = nn_search(G.V',atria,(1:G.nV)',BNN+1,-1,0.0);
% fullPhi = sparse(repmat(1:G.nV,1,BNN+1),idx,exp(-dist.^2/bandwidth),G.nV,G.nV);
[idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
fullPhi = sparse(repmat(1:G.nV,1,BNN+1),idx,exp(-dist.^2/bandwidth),G.nV,G.nV);
fullPhi = (fullPhi+fullPhi')/2;

% PDistMat = squareform(pdist(G.V'));
% fullPhi = exp(-PDistMat.^2/bandwidth);

disp('Constructing full kernel......');
tic;
fullMatProd = fullPhi * sparse(1:G.nV,1:G.nV,Lambda,G.nV,G.nV) * fullPhi;
disp(['full kernel constructed in ' num2str(toc) ' sec.']);

KernelTrace = diag(fullMatProd);
GPLmkIdx = zeros(1,numLmk);

invKn = zeros(numLmk);

cback = 0;
for j=1:numLmk
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d\n',j);
    
    if j == 1
        ptuq = KernelTrace;
    else
        if j == 2
            invKn(1:(j-1),1:(j-1)) = 1/fullMatProd(GPLmkIdx(1),GPLmkIdx(1));
            ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
                .*(invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:)),1)';
        else
            p = fullMatProd(GPLmkIdx(1:(j-2)),GPLmkIdx(j-1));
            mu = 1./(fullMatProd(GPLmkIdx(j-1),GPLmkIdx(j-1))-p'*invKn(1:(j-2),1:(j-2))*p);
            invKn(1:(j-2),1:(j-1)) = invKn(1:(j-2),1:(j-2))*[eye(j-2)+mu*(p*p')*invKn(1:(j-2),1:(j-2)),-mu*p];
            invKn(j-1,1:(j-1)) = [invKn(1:(j-2),j-1)',mu];
            productEntity = invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:);
            ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
                .*productEntity,1)';
        end
    end
    [~,maxUQIdx] = max(ptuq);
    GPLmkIdx(j) = maxUQIdx;
end

p = fullMatProd(GPLmkIdx(1:(end-1)),GPLmkIdx(end));
mu = 1./(fullMatProd(GPLmkIdx(end),GPLmkIdx(end))-p'*invKn(1:(end-1),1:(end-1))*p);
invKn(1:(end-1),:) = invKn(1:(end-1),1:(end-1))*[eye(numLmk-1)+mu*(p*p')*invKn(1:(end-1),1:(end-1)),-mu*p];
invKn(end,:) = [invKn(1:(end-1),end)',mu];
ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx)'...
    .*(invKn*fullMatProd(GPLmkIdx,:)),1)';

end

