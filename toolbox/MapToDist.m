function [Dist, R, T] = MapToDist(V1,V2,map12,VertArea1)
%MAPTODIST Summary of this function goes here
%   Detailed explanation goes here

if size(V1,1)<size(V1,2)
    V1 = V1';
end
if size(V2,1)<size(V2,2)
    V2 = V2';
end

[R,T,Dist] = findBestRigidMotion(V1,V2(map12,:),VertArea1);

end

