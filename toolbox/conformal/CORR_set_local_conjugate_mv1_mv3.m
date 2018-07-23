function [u_star_mv3] = CORR_set_local_conjugate_mv1_mv3(u_star_mv1,v1,v2,v3,u1,u2,u3)

alpha1 = myangle(v2-v1, v3-v1);
%alpha2 = myangle(v3-v2, v1-v2);
alpha3 = myangle(v1-v3, v2-v3);
u_star_mv3 = u_star_mv1 + ...
    0.5*((u1-u2)*cot(alpha3) + (u3-u2)*cot(alpha1) );%the minus because 1->2 is inverse direction





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
%du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );