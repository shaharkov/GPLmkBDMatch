function [R, T, err] = TEETH_find_best_rigid_motion(P,Q,VertArea)
%give two corresponding point sets P,Q Nx3, find best rigid alignment R
%(rotation ) and T (translation), and err the mean-squared-deviation

if size(VertArea,1)<size(VertArea,2)
    VertArea = VertArea';
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
err = sqrt(VertArea'*sum((tP-Q).^2,2));

end