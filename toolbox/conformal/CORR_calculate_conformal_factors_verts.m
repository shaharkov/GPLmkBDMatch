function [LambdaV] = CORR_calculate_conformal_factors_verts(F,Lambda)
% calculate the scaling factors at the verts using average of those at the 
% adjacent triangles

ring = CORR_compute_vertex_face_ring(F);


LambdaV=zeros(length(ring),1);

for k=1:length(ring);
   f = ring{k};
    
   LambdaV(k) = sum(Lambda(f))/length(f);
end



