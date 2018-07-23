function [R, T, err] = findBestRigidMotion(P,Q,w)
%give two corresponding point sets P,Q Nx3, find best rigid alignment R
%(rotation ) and T (translation), and err the mean-squared-deviation

if size(w,1)<size(w,2)
    w = w';
end

n = size(P,1);
if(size(Q,1) ~= n)
    disp('error, P,Q not equal size');
    return;
end
mP = mean(P,1);
mQ = mean(Q,1);

T = mQ-mP;

P = P - ones(n,1)*mP;
Q = Q - ones(n,1)*mQ;

[U,~,V] = svd(P'*Q);
R = V*U';

tP = (R*P')';
err = sqrt(w'*sum((tP-Q).^2,2));

end
