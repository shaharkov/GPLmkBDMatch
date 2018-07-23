function K = LSCM2d(G,bIdx_x,bValues_x,bIdx_y,bValues_y,Lnorm)
V=G.V; %vertices
F=G.F; %faces

if size(V,1)==3 && any(V(3,:))
    error('Graph is not planar');
end
Nv=size(V,2); % number of vertices
Nf=size(F,2); %number of faces
Nb_x=numel(bIdx_x);
Nb_y=numel(bIdx_y);

[gradx,grady] = ComputeGradientMatrix(G);
rot = sparse(1:2*Nf,[Nf+1:2*Nf,1:Nf],[-ones(1,Nf),ones(1,Nf)],2*Nf,2*Nf);
M=[gradx;grady];

A=[M,rot*M];

%% boundary conditions
if ~isempty(bIdx_x) && ~isempty(bIdx_y)
    C_pos_x=sparse(1:Nb_x,bIdx_x,ones(1,Nb_x),Nb_x,Nv);
    C_pos_y=sparse(1:Nb_y,bIdx_y,ones(1,Nb_y),Nb_y,Nv);
elseif ~isempty(bIdx_x) && isempty(bIdx_y)
    C_pos_x=sparse(1:Nb_x,bIdx_x,ones(1,Nb_x),Nb_x,Nv);
    C_pos_y=ones(1,Nv); Nb_y=1; bValues_y=1;
elseif isempty(bIdx_x) && ~isempty(bIdx_y)
    C_pos_x=ones(1,Nv); Nb_x=1; bValues_x=1;
    C_pos_y=sparse(1:Nb_y,bIdx_y,ones(1,Nb_x),Nb_y,Nv);
else
    error('No boundary conditions!');
end
C_pos=[C_pos_x,zeros(Nb_x,Nv);zeros(Nb_y,Nv),C_pos_y];

C=C_pos;

b=[zeros(size(A,2),1);bValues_x(:);bValues_y(:)];

%% L_2 norm
if Lnorm==2
    D=[A'*A,C';C,zeros(size(C,1))];
    x=D\b;

    x=x(1:size(A,2));
end
%% L_inf norm
if strcmp(Lnorm,'inf')
    b=b(2*Nv+1:end);
    cvx_begin
        variable x(2*Nv,1)
        minimize( norm( A * x, Inf ) )
        subject to
            C * x == b
    cvx_end
end

%% result
U=[reshape(x,Nv,2)';zeros(1,Nv)];
K=Mesh('VF',U,F);