function [pV] = CORR_flatten_mesh_by_COT_using_midedge(V,F,M, mV,mF,pmV, cut_face_ind)
%use cot-weights on the mesh to flatten it w.r.t. the flattened mid-edge mesh

[L] = CORR_compute_laplacian_tension(V,F);

%for each vertex calculate its MVC in the space
%and use it to define its planar location
pV = zeros(size(V,1),2);
for k=1:size(V,1)
    v=V(k,:);
    Nind = find(M(k,:));%indices of the k-th vertex's neighbors
    MNind = M(k,Nind);
    %get the coordinate weight from the laplace matrix
    W=-L(k,Nind);
    W = W./L(k,k);    
    pV(k,:)=W*pmV(MNind,:);
end


%for the vertices of the cut_face decide by mid-egde->vertices relation in
%that face

mv = pmV(mF(cut_face_ind,:),:);
A = 0.5*[1 0 1; 1 1 0; 0 1 1];
B = [mv(2,:) ; mv(3,:); mv(1,:)];
v = A\B;
%arrange the order
mf = mF(cut_face_ind,:);
f = F(cut_face_ind,:);
for k=1:3
    if(k==3)
        k1=1;
        k2=2;
    elseif(k==2)
        k1=3;
        k2=1;
    else
        k1=2;
        k2=3;
    end
    ind1 = M(f(k),f(k1));
    ind2 = M(f(k),f(k2));
    find1 = find(mf == ind1);
    find2 = find(mf == ind2);
    m(k) = CORR_find_not_common([1 2 3], [find1 find2]);
end
[sorted per] = sort(m);

pV(F(cut_face_ind,per),:) = v;

