function [A,b] = PlanarThreePtsDeform(X,Y)
%PLANARTHREEPTSDEFORM Summary of this function goes here
%   Detailed explanation goes here

if (size(X,1)>size(X,2))
    X = X';
end
if (size(Y,1)>size(Y,2))
    Y = Y';
end
if (size(X,2)~=3)||(size(Y,2)~=3)
    error('This function only works for three-point linear interpolation!');
end
if (size(X,1)~=2)||(size(Y,1)~=2)
    error('This function only works for planar points!');
end

ProjA = Y/[X;ones(1,size(X,2))];
A = ProjA(:,1:2);
b = ProjA(:,3);


end

