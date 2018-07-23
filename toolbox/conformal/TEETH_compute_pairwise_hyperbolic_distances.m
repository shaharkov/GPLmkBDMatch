function [P] = TEETH_compute_pairwise_hyperbolic_distances(s1,s2)

n1 = length(s1);
n2 = length(s2);

XmY = s1 * ones(1,n2) - ones(n1,1) * s2.';
XY = (s1 * ones(1,n2)) .* (ones(n1,1) * s2') ;
 
% Hyperbolic distance
P = atanh( abs(XmY)./abs(ones(size(XY)) - XY ) );
