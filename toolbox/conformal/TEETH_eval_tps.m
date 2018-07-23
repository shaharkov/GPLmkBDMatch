function [fP] = TEETH_eval_tps(tpsf,P)

n=size(P,1);
cen = tpsf.centers';
m = size(cen,1);

%create the collocation matrix CM
PX = P(:,1)*ones(1,m)-ones(n,1)*cen(:,1)';
PY = P(:,2)*ones(1,m)-ones(n,1)*cen(:,2)';

CM = PX.^2+PY.^2;
CM = CM.*log(CM);
%tind = 1:(size(CM,1)+1): prod(size(CM));
CM(isnan(CM))=0;

%evaluate
fP = CM*tpsf.coefs(1,1:m)';
%linear part
fP = fP + P(:,1)*tpsf.coefs(1,m+1) + P(:,2)*tpsf.coefs(1,m+2) + ...
    ones(n,1)*tpsf.coefs(1,m+3);

for k=2:size(tpsf.coefs,1)
    fP(:,k) = CM*tpsf.coefs(k,1:m)';
    %linear part
    fP(:,k) = fP(:,k) + P(:,1)*tpsf.coefs(k,m+1) + P(:,2)*tpsf.coefs(k,m+2) + ...
        ones(n,1)*tpsf.coefs(k,m+3);
end
% 
% if (size(tpsf.coefs,1)>1)
%     fP(:,2) = CM*tpsf.coefs(2,1:m)';
%     %linear part
%     fP(:,2) = fP(:,2) + P(:,1)*tpsf.coefs(2,m+1) + P(:,2)*tpsf.coefs(2,m+2) + ...
%         ones(n,1)*tpsf.coefs(2,m+3);
% end
% 
% if (size(tpsf.coefs,1)>2)
%     fP(:,3) = CM*tpsf.coefs(3,1:m)';
%     %linear part
%     fP(:,3) = fP(:,3) + P(:,1)*tpsf.coefs(3,m+1) + P(:,2)*tpsf.coefs(3,m+2) + ...
%         ones(n,1)*tpsf.coefs(3,m+3);
% end
% 



