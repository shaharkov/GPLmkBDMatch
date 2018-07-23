function [K,Kav] = SNCCM2d(G,bIdx_x,bValues_x,bIdx_y,bValues_y,Lnorm,midedge)
% Symmetric Non-Conforming Conformal Map

if size(G.V,1)==3 && any(G.V(3,:))
    error('Graph is not planar');
end

Nf=size(G.F,2);
Nv=size(G.V,2);
Nb_x=numel(bIdx_x);
Nb_y=numel(bIdx_y);

[K,V2V,E2E]=G.Soupalize;
[gradx,grady] = K.ComputeGradientMatrix;
M=[gradx;grady];


%constrain othogonality
rot = sparse(1:2*Nf,[Nf+1:2*Nf,1:Nf],[-ones(1,Nf),ones(1,Nf)],2*Nf,2*Nf); % rotation matrix for each face
C_orth=[M,rot*M];


%constrain boundary
if ~isempty(bIdx_x) && ~isempty(bIdx_y)
    C_pos_x=sparse(1:Nb_x,1:Nb_x,1./sum(V2V(bIdx_x,:),2),Nb_x,Nb_x)*V2V(bIdx_x,:); 
    C_pos_y=sparse(1:Nb_y,1:Nb_y,1./sum(V2V(bIdx_y,:),2),Nb_y,Nb_y)*V2V(bIdx_y,:); 
elseif ~isempty(bIdx_x) && isempty(bIdx_y)
    C_pos_x=sparse(1:Nb_x,1:Nb_x,1./sum(V2V(bIdx_x,:),2),Nb_x,Nb_x)*V2V(bIdx_x,:); 
    C_pos_y=ones(1,Nv)*V2V; Nb_y=1; bValues_y=1;
elseif isempty(bIdx_x) && ~isempty(bIdx_y)
    C_pos_x=ones(1,Nv)*V2V; Nb_x=1; bValues_x=1;
    C_pos_y=sparse(1:Nb_y,1:Nb_y,1./sum(V2V(bIdx_y,:),2),Nb_y,Nb_y)*V2V(bIdx_y,:); 
else
    error('No boundary conditions!');
end
C_pos=[C_pos_x,zeros(Nb_x,3*Nf);zeros(Nb_y,3*Nf),C_pos_y];

% constrain mid-edges
if exist('midedge','var')==1 && midedge==1
    Ic_K=K.ComputeV2E;
    C_mid=E2E(sum(E2E,2)==0,:)*abs(Ic_K/2); C_mid=kron(eye(2),C_mid);
else
    C_mid=[];
end
C=[C_orth;C_mid;C_pos];

b=[zeros(6*Nf,1);zeros(size(C_orth,1),1);zeros(size(C_mid,1),1);bValues_x(:);bValues_y(:)];

A=speye(3*Nf)-V2V'*sparse(1:Nv,1:Nv,1./sum(V2V,2),Nv,Nv)*V2V;
%A(bIdx,:)=[];
A=blkdiag(A,A);

%% L_2 norm
if Lnorm==2
%    x=pinv(full([A'*A,C';C,zeros(size(C,1))]))*b;
    D=[A'*A,C';C,zeros(size(C,1))];
    x=D\b;

    x=x(1:size(A,2));
end
%% L_inf norm
if strcmp(Lnorm,'inf')
    b=b(6*Nf+1:end);
    cvx_begin
        variable x(6*Nf,1)
        minimize( norm( A * x, Inf ) )
        subject to
            C * x == b
    cvx_end
end

%% Create new mesh
K.V(1:2,:)=reshape(x,[],2)'; K.V(3,:)=0;
K=Mesh('VF',K.V,K.F);
if nargout>1
    Kav=Mesh(G);
    Kav.V(1:2,:)=(sparse(1:Nv,1:Nv,1./sum(V2V,2),Nv,Nv)*V2V*K.V(1:2,:)')';
    Kav.V(3,:)=0;
end