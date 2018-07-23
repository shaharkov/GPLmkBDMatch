function [pmV ] = CORR_map_mesh_to_plane_nonconforming(V,F,mF,seed_face,M,E2V,numE, reflect_mesh)
% map the mesh V,F to plane using the non-conforming method of polthier
% use the seed face to cut the mesh open
% output: the midpoint mesh and its embedding to plane

oF = F; %rememeber the face list for later

%cut out the seed face
ioutf = seed_face;
outf = F(ioutf,:);
outf=sort(outf);
F(ioutf,:)=[];

[L] = CORR_compute_laplacian_tension(V,F);
L1 = L;
outf = sort(outf);
L1(outf(1),:) = [];
L1(outf(2)-1,:) = [];

L1rows = size(L1,1);

L1(L1rows+1,outf(1)) = 1;
L1(L1rows+2,outf(2)) = 1;

b = zeros(L1rows,1);
b(L1rows+1) = -1;
b(L1rows+2) = 1;
u = L1\b;

% %with linear system
% [e_u_star] = CORR_calculate_conjugate_harmonic(F,V,u,M,E2V,numE);
%withOUT linear system
imissing_f = seed_face;
[e_u_star] = CORR_calculate_conjugate_harmonic_faster(oF,V,mF,u,M,E2V,numE,imissing_f);

[mu] = CORR_get_midpoint_values_of_function(oF,u, M, E2V, numE);

if(reflect_mesh==0)
    pmV = [mu e_u_star]; 
else
    pmV = [mu -e_u_star]; 
end

