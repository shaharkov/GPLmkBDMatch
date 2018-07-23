function [dpmV] = CORR_transform_to_disk_new_with_dijkstra(pmV,mF,e2v,center_ind)
%D is the graph dijkstra distance - use for deciding on vertices near branching
%points

[mBoundary, ~] = CORR_locate_midedge_boundary_vertices(mF);

%we want to take only the largest boundary
[bcc] = CORR_find_midedge_boudnary_connected_component(mF,e2v);
cc = CORR_sort_and_delete_repeats(bcc);
cc_count= zeros(1,length(cc));
for j=1:length(cc)
    cc_count(j) = length(find(bcc == cc(j)));    
end
[~, maxind] = max(cc_count);
mBoundary = mBoundary(bcc==cc(maxind));

pmB = pmV(mBoundary,:);

[Xmin, imin] = min(pmB(:,1));
[Xmax, imax] = max(pmB(:,1));
cen = 0.5*(pmB(imin,:)+pmB(imax,:));
%move to center X axis and normalize
pmV = (pmV - ones(size(pmV,1),1)*cen)*4/(Xmax-Xmin);
% pmB = pmV(mBoundary,:);

%transform to unit disk all except the boundary
P = pmV(:,1) + 1i*pmV(:,2);
temp = 0.5.*(P - sqrt(P.^2 - 4));
onebranch = temp;
secondbranch = 0.5.*(P + sqrt(P.^2 - 4));
temp = abs(temp);

%nan all points with temp close to 1
eps = 0.01;  %was 0.01 and worked ok. consider to change back if fails often
tind = find(abs(temp-1)<eps);
P(tind)=NaN;
temp(tind)=NaN; %for later in the loop below we use it to find the closest no-nan to nan points
P = (temp <= 1).*0.5.*(P - sqrt(P.^2 - 4)) + (temp > 1).*0.5.*(P + sqrt(P.^2 - 4));

%%now use the points transformed to set points undecided before (temp<eps)
A = triangulation2adjacency(mF);
[D] = TEETH_compute_distance_graph(A,tind);


%for each tind find tes closest decided neighbor and set it with the same
%sign of temp-1
for j=1:size(D,1)
    progressbar(j,size(D,1),40);
    td = D(j,:);
    [sa, si] = sort(td);
    stemp = P(si);
    stemp=TEETH_remove_nans(stemp);
    %find the correct branch to use for the points near the branching
    if (abs(onebranch(tind(j))-stemp(1)) < abs(secondbranch(tind(j))-stemp(1)) )
        P(tind(j)) = onebranch(tind(j));
    else
        P(tind(j)) = secondbranch(tind(j));
    end
end


% %further take the center point to the center by disc mobius
center=P(center_ind);
mob = [ 1 -center ; -conj(center) 1 ];
[P] = CORR_apply_moebius_as_matrix(mob,P);

dpmV = [real(P) imag(P)];

