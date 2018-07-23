function [Lambda] = CORR_calculate_conformal_factors(F,V,pV)
% calculate the triangles scaling in the flattening v->pV
%that is: lambda(f) = area(f)/area(pf), where pf denotes the flattened face


Lambda=zeros(size(F,1),1);

for k=1:size(F,1);
   f = F(k,:);
   
   v1 = V(f(1),:);
   v2 = V(f(2),:);
   v3 = V(f(3),:);
   A = [v2-v1;v3-v1];
   
   u1 = pV(f(1),:);
   u2 = pV(f(2),:);
   u3 = pV(f(3),:);
   B = [u2-u1;u3-u1];
   
   Lambda(k) = sqrt(abs(det(A*A'))/abs(det(B*B')));
end



