function [bcc] = CORR_find_midedge_boudnary_connected_component(mF1,E2V_1)
%find for each boundary *mid-edge* vertex a connected component
%the output is indices with indices of connected components

[mBoundary ring] = CORR_locate_midedge_boundary_vertices(mF1);
 
bcc = 1:length(mBoundary); %each boundary mvertex is in its own connected component

%all mid-edge boundary vertices tounch two vertices
mev_2_v = E2V_1(mBoundary,:);

%for every boundary me-v : make its component the same as all other who
%share a vertex with him
for k=length(mBoundary):-1:1
    v = mev_2_v(k,:);
    tind = find(mev_2_v(:,1) == v(1) | mev_2_v(:,1) == v(2) | ...
        mev_2_v(:,2) == v(1) | mev_2_v(:,2) == v(2));
    aind=[];
    for j=1:length(tind)
        aind = [aind find(bcc == bcc(tind(j))) ];
    end
    bcc(aind) = k;
end


