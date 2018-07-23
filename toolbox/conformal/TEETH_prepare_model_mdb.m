function [Mdb1] = TEETH_prepare_model_mdb(model_filename, target_filename)
% This function takes as input an .off mesh and outpus mdb file
% that is used later for comparing pairs of surfaces.

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp(' ');
disp('>>----------------------------------------------------------------');
disp('Preparing mesh:');
disp(model_filename);



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%parameters (independent) file -> to set parameters open this file
% parameters;
parameters_gui;

%add paths to the filenames
% model_filename = [meshes_path model_filename];
% target_filename =[mdbs_path target_filename];


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%read model (now only .off files)
[V1,F1]=read_off( model_filename );
Mdb1=options;
Mdb1.V1=V1'; Mdb1.F1=F1';

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Normalize the surface area to 1
area1 =  CORR_calculate_area(Mdb1.F1,Mdb1.V1);
Mdb1.V1 = Mdb1.V1*sqrt(1/area1);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PREPARE MESH: SPREAD POINTS AND FLATTEN
% %spread based on geodesics
% [Mdb1] =  CORR_prepare_and_flatten(Mdb1);
%spread based on Euclidean distances (good stable approximation)
[Mdb1] =  CORR_prepare_and_flatten_EUC(Mdb1);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% map domain to disk (add the infinity point back as sample point)
% transfer the indices of center point to the mid-edge
center_ind = Mdb1.center_ver_1;
tind = Mdb1.V2F_1{center_ind};%v_max_V1(kk);
%faces have the same index in mid-edge and original mesh
center_ind = Mdb1.mF1(tind(1),1);

%%older but faster:
%  [Mdb1.dpmV1] = CORR_transform_to_disk_new(Mdb1.pmV1,Mdb1.mF1, Mdb1.E2V_1, center_ind);
% % % %new: use dijkstra to map domain to disk (we assume disk-type surfaces)
[Mdb1.dpmV1] = CORR_transform_to_disk_new_with_dijkstra( ...
   Mdb1.pmV1,Mdb1.mF1, Mdb1.E2V_1, center_ind);

%%add the missing face
Mdb1.pmF1 = Mdb1.mF1;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Map the original mesh to the disk using the mid-edge structure
disp('Flattening the ORIGINAL mesh using the mid-edge flattening...')
Mdb1.pV1 = CORR_flatten_mesh_by_COT_using_midedge(Mdb1.V1,Mdb1.F1,...
    Mdb1.M1, Mdb1.mV1, Mdb1.mF1, Mdb1.dpmV1, Mdb1.cut_face_1);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% remove dulicate points and indices of NaN points from density points (due to second conformal
% map). delete NaN point, if appeared
dind = find(isnan(Mdb1.dpmV1(Mdb1.me_density_points,1)));
if( ~isempty(dind) )
    disp(['--> found ' num2str(length(dind)) ' NaN sample points, deleting them.'])
    Mdb1.me_density_points(dind)=[];
%     Mdb1.samples_GeoField(dind,:)=[];
    Mdb1.num_of_spread_points = Mdb1.num_of_spread_points - length(dind);
end

% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % delete boundary points
% geo = Mdb1.samples_GeoField;
% dist_to_boundary = min(geo(:,Mdb1.vBoundary),[],2);
% Mdb1.diameter = max(max(geo));
% del_ind = find(dist_to_boundary < minimal_distance_to_boundary * Mdb1.diameter);
% Mdb1.num_of_spread_points = Mdb1.num_of_spread_points - length(del_ind);
% Mdb1.samples_GeoField(del_ind,:) = [];
% Mdb1.sample_pnts1(del_ind)=[];
% % Mdb1.samples_geodesic_mat(del_ind,:)=[];
% % Mdb1.samples_geodesic_mat(:,del_ind)=[];

% % %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % % take sample points far away from the boundary
% % geo = Mdb1.samples_GeoField;
% % dist_to_boundary = min(geo(:,Mdb1.vBoundary),[],2);
% % Mdb1.diameter = max(max(geo));
% % del_ind = find(dist_to_boundary < minimial_distance_to_boundary_of_compared * Mdb1.diameter);
% % Mdb1.boundary_sample_pnts1 = del_ind;% indices in sample points. %old: Mdb1.sample_pnts1(del_ind);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Calculate the conformal factors
% first, for the me-faces and then to the me-vertices.
[Mdb1.Lambdam] = CORR_calculate_conformal_factors(Mdb1.pmF1,Mdb1.mV1,Mdb1.dpmV1);
[Mdb1.LambdamV] = CORR_calculate_conformal_factors_verts(Mdb1.pmF1,Mdb1.Lambdam);

% % %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % % CLAMP the conformal factors
% disp('Clamp high conformal factors...!')
% tind = Mdb1.LambdamV> 5;
% Mdb1.LambdamV(tind)=5;
% tind = Mdb1.Lambdam> 5;
% Mdb1.Lambdam(tind)=5;
% %     for kk=1:3
% %         [mE1] = CORR_average_face_data_over_vertex_ring(Mdb1.LambdamV,Mdb1.mF1);
% %         Mdb1.LambdamV = CORR_average_over_ring(Mdb1.mV1, Mdb1.mF1, mE1)';
% %     end

% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Choose generating set (feature point) for checking transformations  
disp('Building generating set...');
switch Mdb1.method_of_feature_pick
    case 'extrema_conf_factors'
        %%compute local maxima of conformal factors.
        Mdb1.mind1 = TEETH_find_local_maxima_of_conformal_factors(Mdb1);
        
    case 'gauss_curvature'        
        type=1;
        
    case 'gauss_maxima_curvature'
        type=1.1;
            
    case 'mean_curvature'        
        type=2;
        
    case 'gauss_and_mean_curvature'        
        type=3;
end

if(~strcmp(Mdb1.method_of_feature_pick,'extrema_conf_factors'))
switch type
    case {1,1.1,2,3}
        mind1 = TEETH_find_local_maxima_of_gauss_mean_curvature(Mdb1,type,Mdb1.smooth_fields, Mdb1.local_max_width);

        %move to mid-edge indices
        v_max1 = zeros(length(mind1),1);
        for kk=1:length(mind1)
            tind = Mdb1.V2F_1{mind1(kk)};%v_max_V1(kk);
            %faces have the same index in mid-edge and original mesh
            v_max1(kk) = Mdb1.mF1(tind(1),1);
        end
        
        Mdb1.mind1 = v_max1;
end
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %remove sample points near the boundary 
disp('Removing generating points near boundary...');
[vBnd] = CORR_locate_boundary_vertices(Mdb1.E2V_1 , Mdb1.V2F_1);
nb = length(vBnd);
nv = length(Mdb1.mind1);
dind=zeros(nv,1);
for jj=1:nb
    dist_to_bnd = Mdb1.mV1(Mdb1.mind1,:) - ones(nv,1)*Mdb1.V1(vBnd(jj),:);
    dist_to_bnd = sum(dist_to_bnd.^2,2);
    dind = dind + (dist_to_bnd < (Mdb1.surfdiam*Mdb1.remove_boundary_portion)^2  );
end
Mdb1.mind1(dind>0)=[];
%

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% find NaN factors (due to transform to disk) and average neighbors
nind = find(isnan(Mdb1.LambdamV)); %all NaN vertices
if(~isempty(nind))
    PP = Mdb1.mV1; %look for closest me-vert on original mesh which is not NaN later
    PP(nind,:)=[]; %remove the ones with NaN
    LL = Mdb1.LambdamV;
    LL(nind)=[];
    tree = KDTreeSearcher(PP);
    idxs = tree.knnsearch(Mdb1.mV1(nind,:));
%     tree = kdtree_build(PP);
%     idxs = kdtree_nearest_neighbor( tree, Mdb1.mV1(nind,:) ); %these are the ones to take the conformal factors from
    Mdb1.LambdamV(nind) = LL(idxs);
end

% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %create cell on parameter domain of first surface as voronoi cells
% %and calculate the area in parameter domain of the cells
% %add the points on the unit disc to get finite voronoi cells
% spind1 = Mdb1.sample_pnts1;
% X = [Mdb1.dpmV1(spind1,1) Mdb1.dpmV1(spind1,2)];
% tet=0:0.05:2*pi; tet=tet';
% X_on_circ = [real(1.2*exp(i*tet)) imag(1.2*exp(i*tet))];
% [VO1 C1] = voronoin([X ; X_on_circ]);
% % voronoi([X(:,1) ; X_on_circ(:,1)],[X(:,2) ; X_on_circ(:,2)]);
% C1 = C1(1:size(X,1));
% %calculate the area of the cells
% areaC1 = zeros(length(C1),1);
% for k=1:length(areaC1)
%     XX = VO1(C1{k}, 1);
%     YY = VO1(C1{k}, 2);
%     areaC1(k) = polyarea(XX,YY);
% end
% Mdb1.areaC1 = areaC1;

% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % transform the conformal factors to hyperbolic
% spind1 = Mdb1.sample_pnts1;
% X1 = [Mdb1.dpmV1(spind1,1) Mdb1.dpmV1(spind1,2)]; %surface 1
% Z=X1(:,1)+1i*X1(:,2);
% mu = Mdb1.LambdamV(spind1);
% Mdb1.LambdamV_samples_hyperbolic  = mu .* (1 - Z.*conj(Z)).^2; %transforming the density to hyperbolic density


% % %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % prepare TPS interpolant of the conformal factors (both w.r.t hyperbolic
% % and Euclidean metrics)
% [Mdb1.tps_mu] = TEETH_build_TPS_for_conformal_factor(Mdb1,options.tps_K,options.tps_p,1);
% Mdb1.tps_nu = Mdb1.tps_mu; %the same thing just depends if this surface is 1st or 2nd
% [Mdb1.tps_mu_EUC] = TEETH_build_TPS_for_conformal_factor(Mdb1,options.tps_K,options.tps_p,2);
% Mdb1.tps_nu_EUC = Mdb1.tps_mu_EUC; %the same thing just depends if this surface is 1st or 2nd


% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %build kd-tree of the sample points
% Mdb1.tree_sam = kdtree_build( [Mdb1.dpmV1(spind1,1) Mdb1.dpmV1(spind1,2)] );
% Mdb1.tree_all = kdtree_build( [Mdb1.dpmV1(:,1) Mdb1.dpmV1(:,2)] );

% % %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % remove unnecessary and large data
% Mdb1.samples_GeoField = []; %%dont save this - too much space!
% 

%%%%%%%%%%%%%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% TG adds: to facilitate display in GUI, rotate each tooth to a "standard
% position"
Mdb1.center_of_mass = mean(Mdb1.V1);
V = Mdb1.V1 - repmat(mean(Mdb1.V1), size(Mdb1.V1, 1), 1);
bdV = CORR_locate_boundary_vertices(Mdb1.E2V_1, Mdb1.V2F_1);

bdV_coord = Mdb1.V1(bdV, :);
bdV_coord_centered = bdV_coord - repmat(mean(Mdb1.V1), size(bdV_coord, 1), 1);
bdV_averaged_normal = mean(-bdV_coord_centered./repmat(sqrt(sum(bdV_coord_centered.^2, 2)), 1, 3));
bdV_averaged_normal = bdV_averaged_normal./repmat(sqrt(sum(bdV_averaged_normal.^2)), 1, 3);
if (det([null(bdV_averaged_normal), bdV_averaged_normal'])>0)
    Mdb1.InvRotateMatrix = [null(bdV_averaged_normal), bdV_averaged_normal'];
else
    temp = null(bdV_averaged_normal);
    temp(:, 1) = -temp(:, 1);
    Mdb1.InvRotateMatrix = [temp, bdV_averaged_normal'];
end
Mdb1.V1_for_display = (Mdb1.InvRotateMatrix\V')';

mV = Mdb1.mV1 - repmat(mean(Mdb1.V1), size(Mdb1.mV1, 1), 1);
Mdb1.mV1_for_display = (Mdb1.InvRotateMatrix\mV')';
[Mdb1.vertex_normal_for_display, Mdb1.face_normal_for_display] = compute_normal(Mdb1.mV1_for_display, Mdb1.mF1);

% %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% save mdb
save(target_filename,'Mdb1');

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
disp('Finished preparing:');
disp(target_filename);
disp('<<----------------------------------------------------------------');



%
% %%%%%
% % plot the mesh and the sigma set and the cut face
% figure
% h1=trimesh(Mdb1.F1,Mdb1.V1(:,1),Mdb1.V1(:,2),...
%     Mdb1.V1(:,3),'EdgeColor','black','FaceColor',[0.8 1 0.8]);
% %set(Mdb.h1,'EdgeAlpha',0.2);
% set(h1,'Edgecolor','none','FaceAlpha',1);
% axis equal
% axis off
% lighting flat
% camlight('headlight');
% %cameramenu(Mdb.fig1);
% hold on
%
% myCol = [1 0 0 ; 0 1 0; 0 0 1; 1 1 0; 1 0 1 ; 0 1 1; 1 1 1];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];myCol = [myCol ; 0.75*myCol];
%
% mV1=Mdb1.mV1;
% %V1=Mdb1.V1;
% scatter3(mV1(Mdb1.sample_pnts1,1),mV1(Mdb1.sample_pnts1,2),mV1(Mdb1.sample_pnts1,3),40, [0 0 1],'filled');
% % scatter3(V1(v_max_V1,1),V1(v_max_V1,2),V1(v_max_V1,3),40, [0 0 1],'filled');
% colorbar
% tempind = Mdb1.cut_ver_1;
% V1=Mdb1.V1;
% scatter3(V1(tempind,1),V1(tempind,2),V1(tempind,3),100, [0 0 0],'filled');
%
% tempind = Mdb1.center_ver_1;
% V1=Mdb1.V1;
% scatter3(V1(tempind,1),V1(tempind,2),V1(tempind,3),100, [1 0 0],'filled');
%



