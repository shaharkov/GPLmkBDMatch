function A = ComputeTriangleAngles(G,EL)

if exist('EL','var')==0
    EL=G.ComputeEdgeLengths;
end

L1=EL(:,1);
L2=EL(:,2);
L3=EL(:,3);
A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
A = [A1,A2,A3];
A = acos(A);

end

