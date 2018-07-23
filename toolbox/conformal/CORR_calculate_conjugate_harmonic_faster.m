function [e_u_star, M, e2v, nume] = CORR_calculate_conjugate_harmonic_faster(F,V,mF,u,M,e2v,nume,imissing_f)
%a procedure to calculate the conjugate harmonic function u_star to u on
%the mid-edges according to Polthier. This is to create conformal mapping
%of the mid edges.
% we use here th dual cot-formula

%construct the mid-edge vertex - mid-edge face ring (which has the same
%face indices as th original mesh
ring = compute_vertex_face_ring(mF);

%traverse the mesh starting for the first face
tobe_visited = zeros(nume,1);
tobe_visited(1) = 1;
tobe_len = 1;
e_u_star = NaN*ones(nume,1);
e_u_star(M(F(tobe_visited(1),1),F(tobe_visited(1),2))) = 0;%set the addditive constant

while(tobe_len>0 ) %while we haven't finished traversing the mesh
    indf = tobe_visited(tobe_len);
    f = F(indf,:);
    tobe_len = tobe_len-1;
    
    %the three mid-edge vertices in this triangle
    imv1 = M(f(2),f(3));
    imv2 = M(f(3),f(1));
    imv3 = M(f(1),f(2));
    
    %         %for debug
    %     if( imv1 == 4 || imv2 == 4 || imv3 == 4)
    %         disp('got to 4');
    %     end
    
    %add to the tobe_visited only the neighboring faces which in the common edge there is no
    %data
    if (isnan(e_u_star(imv1))) %if the mid-edge vertex 1 has no data
        neigh_f=ring{imv1};
        if(neigh_f(1) == indf)
            neigh_f(1)=[];
        end
        if(~isempty(neigh_f))
            if (neigh_f(1) ~= imissing_f)%make sure we are not passing through the missing face
                tobe_len = tobe_len+1;
                tobe_visited(tobe_len) = neigh_f(1); %add this face to be visited
            end
        end
    end
    if (isnan(e_u_star(imv2))) %if the mid-edge vertex 2 has no data
        neigh_f=ring{imv2};
        if(neigh_f(1) == indf)
            neigh_f(1)=[];
        end
        if(~isempty(neigh_f))
            if (neigh_f(1) ~= imissing_f)
                tobe_len = tobe_len+1;
                tobe_visited(tobe_len) = neigh_f(1); %add this face to be visited
            end
        end
    end
    if (isnan(e_u_star(imv3))) %if the mid-edge vertex 3 has no data
        neigh_f=ring{imv3};
        if(neigh_f(1) == indf)
            neigh_f(1)=[];
        end
        if(~isempty(neigh_f))
            if (neigh_f(1) ~= imissing_f)
                tobe_len = tobe_len+1;
                tobe_visited(tobe_len) = neigh_f(1); %add this face to be visited
            end
        end
    end
    
    if (~isnan(e_u_star(imv1))) %if the mid-edge vertex 1 has data
        %dont change f
    elseif (~isnan(e_u_star(imv2)))
        f = [f(2) f(3) f(1)];
    elseif (~isnan(e_u_star(imv3)))
        f = [f(3) f(1) f(2)];
    else
        %not error - can happen
        %         disp('Error 2341 - stopping...')
        %         return
    end
    imv1 = M(f(2),f(3)); %set the midedges vertices indices again
    imv2 = M(f(3),f(1));
    imv3 = M(f(1),f(2));
    
    %     %for debug
    %     e_u_star(imv3) = 1;
    %     e_u_star(imv2) = 1;
    %     e_u_star(imv1) = 1;
    %     if(tobe_len == 0)
    %         disp('for debug');
    %         break;
    %     end
    
    
    v1 = V(f(1),:); %the regular vertices
    v2 = V(f(2),:);
    v3 = V(f(3),:);
    u1 = u(f(1)); %the discrete harmonic values st the vertices v_i
    u2 = u(f(2)); %the discrete harmonic values
    u3 = u(f(3)); %the discrete harmonic values
    
    %set the u_star at the other two vertices
    %         if(reflect_mesh == 0)
    e_u_star(imv3) = CORR_set_local_conjugate_mv1_mv3(e_u_star(imv1),v1,v2,v3,u1,u2,u3);
    e_u_star(imv2) = CORR_set_local_conjugate_mv1_mv3(e_u_star(imv3),v3,v1,v2,u3,u1,u2);
    %         else %should reflect the flattening
    %             e_u_star(imv3) = -CORR_set_local_conjugate_mv1_mv3(-e_u_star(imv1),v1,v2,v3,u1,u2,u3);
    %             e_u_star(imv2) = -CORR_set_local_conjugate_mv1_mv3(-e_u_star(imv3),v3,v1,v2,u3,u1,u2);
    %         end
    
    %     %for debug
    %     if( imv2 == 4 || imv3 == 4)
    %         disp('got to 4');
    %     end
    
end
