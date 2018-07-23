function [K, Kav] = NCCM2d(G,bIdx,bValues)
% Non-Conforming Conformal Map (a-la Pinkal-Polthier)
if size(G.V,1)==3 && any(G.V(3,:))
    error('Graph is not planar');
end

Nf=size(G.F,2);
Nv=size(G.V,2);
Nb=numel(bIdx);

[K,V2V,E2E]=G.Soupalize;
[gradx,grady] = K.ComputeGradientMatrix;
M=[gradx;grady];

% find discrete harmonic part
L=G.ComputeLaplacian;
C_pos=sparse(1:Nb,bIdx,ones(1,Nb),Nb,Nv); 
Lr=L(setdiff(1:Nv,bIdx),:);
u=[Lr'*Lr,C_pos';C_pos,zeros(size(C_pos,1))]\[zeros(Nv,1);bValues(:)];
u=u(1:Nv);


% find conjugate part
%constrain othogonality
rot = sparse(1:2*Nf,[Nf+1:2*Nf,1:Nf],[-ones(1,Nf),ones(1,Nf)],2*Nf,2*Nf); % rotation matrix for each face
C_orth=M;
b_orth=-rot*M*V2V'*u;

% constrain mid-edges
Ic_K=K.ComputeV2E;
C_mid=E2E(sum(E2E,2)==0,:)*abs(Ic_K/2); 

C=[C_orth;C_mid;ones(1,3*Nf)];

b=[b_orth;zeros(size(C_mid,1),1);0];

v=C\b;
v=v(1:3*Nf);

%% Create new mesh
K=Mesh('VF',[V2V'*u,v,zeros(3*Nf,1)]',K.F);
if nargout>1
    Kav=Mesh(G);
    Kav.V(1:2,:)=(sparse(1:Nv,1:Nv,1./sum(V2V,2),Nv,Nv)*V2V*K.V(1:2,:)')';
    Kav.V(3,:)=0;
end