function K = NCCM(G,bIdx,bValues,Lnorm)
% Non-Conforming Conformal Map (a-la Pinkal-Polthier)
[K,V2V]=G.Soupalize;
Ic_G = G.ComputeV2E;
Ic_K = K.ComputeV2E;
% This find the parent j of edge i in the soup
E2E=mod(find(abs(Ic_G'*V2V*Ic_K)==2),size(Ic_G,2)); E2E(E2E==0)=size(Ic_G,2);
E2E=sparse(1:size(Ic_K,2),E2E,1,size(Ic_K,2),size(Ic_G,2)); E2E=E2E';
% index of last non-zero element of each row of E2E
j = sum(cumsum((E2E(:,end:-1:1) ~= 0), 2) ~= 0, 2);
E2E(sub2ind(size(E2E),1:size(E2E,1),j'))=-1;
Nf=size(G.F,2);
Nv=size(G.V,2);

%Fn=G.Fn; % G.Fn==K.Fn
Fn=K.ComputeFaceNormals;
[gradx,grady,gradz] = K.ComputeGradientMatrix;
M=sparse([gradx;grady;gradz]);

% sparsify this part
rot=zeros(3*Nf);
i=1;
Fni=Fn(:,i); % current face normal
rot(i+[0,Nf,2*Nf],i+[0,Nf,2*Nf])=[0,-Fni(3), Fni(2);Fni(3),0,-Fni(1);-Fni(2),Fni(1),0];
for i=2:Nf
    Fni=Fn(:,i); % current face normal
    rot(i+[0,Nf,2*Nf],i+[0,Nf,2*Nf])=[0,-Fni(3), Fni(2);Fni(3),0,-Fni(1);-Fni(2),Fni(1),0];
end
rot=sparse(rot);
% end sparsify

C_orth=[M,-rot*M];
%C_orth(
C_pos=sparse(1:2*numel(bIdx),[bIdx;bIdx+3*Nf],1,2*numel(bIdx),6*Nf);
C_mid=E2E(sum(E2E,2)==0,:)*abs(Ic_K/2); C_mid=kron(eye(2),C_mid);
C=[C_orth;C_mid;C_pos];

b=[zeros(6*Nf,1);zeros(3*Nf,1);zeros(size(C_mid,1),1);reshape(bValues',[],1)];

A=speye(3*Nf)-V2V'*sparse(1:Nv,1:Nv,1./sum(V2V,2),Nv,Nv)*V2V;
A(bIdx,:)=[];
A=blkdiag(A,A);
%% L_2 norm
if Lnorm==2
%    x=pinv(full([A'*A,C';C,zeros(size(C,1))]))*b;
    x=[A'*A,C';C,zeros(size(C,1))]\b;

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


K.V(1:2,:)=reshape(x,[],2)'; K.V(3,:)=0;
K=MeshGraph('VF',K.V,K.F);