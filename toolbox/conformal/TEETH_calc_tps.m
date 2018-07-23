function [ftps] = TEETH_calc_tps(cen,F)
%calculate the coefficients of the tps interpolant phi(P_i)=F_i

m = size(cen,1);

%create the collocation matrix CM
PX = cen(:,1)*ones(1,m)-ones(m,1)*cen(:,1)';
PY = cen(:,2)*ones(1,m)-ones(m,1)*cen(:,2)';

CM = PX.^2+PY.^2;
CM = CM.*log(CM);
CM(isnan(CM))=0;

L = [cen ones(m,1)];

M = [CM L; L' zeros(3,3)];

siz = size(F);
%for debug
if(rcond(M)>10^6)
    ftps=0;
    return;
end
coefs = M\[F ; zeros(3,siz(2))];


ftps.centers = cen';
ftps.coefs = coefs';

