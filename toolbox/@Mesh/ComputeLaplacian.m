function [L,M]=ComputeLaplacian(G)

V=G.V';
F=G.F';
Nv=size(V,1);
L1 = sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2));
L2 = sqrt(sum((V(F(:,1),:)-V(F(:,3),:)).^2,2));
L3 = sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2));
S=(L1+L2+L3)/2;
Ar=sqrt(S.*(S-L1).*(S-L2).*(S-L3));
A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
A = [A1,A2,A3];
A = acos(A);

I = [F(:,1);F(:,2);F(:,3)];
J = [F(:,2);F(:,3);F(:,1)];
S = 0.5*cot([A(:,3);A(:,1);A(:,2)]);
In = [I;J;I;J];
Jn = [J;I;I;J];
Sn = [-S;-S;S;S];
L = sparse(In,Jn,Sn,Nv,Nv);

M=sparse(1:Nv,1:Nv,1./(G.F2V'*Ar/3),Nv,Nv);