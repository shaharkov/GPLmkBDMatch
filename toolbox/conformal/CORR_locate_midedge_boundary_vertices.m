function [mBoundary ring] = CORR_locate_midedge_boundary_vertices(mF)

ring = CORR_compute_vertex_face_ring(mF);

nver = max(max(mF));

nmax = floor(nver/3);
mBoundary =zeros(1,nmax);
end_ind=0;
for k=1:nver
    if (length(ring{k}) < 2)
        end_ind = end_ind+1;
        mBoundary(end_ind) = k;
    end
    
    
end

mBoundary(end_ind+1:nmax) =[];