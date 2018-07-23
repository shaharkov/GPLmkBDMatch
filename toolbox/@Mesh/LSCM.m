function K = LSCM(G,bIdx,bValues)
V=G.V; %vertices
F=G.F; %faces
Nv=size(V,2); % number of vertices
Nf=size(F,2); %number of faces
Fn=G.ComputeFaceNormals; %face normals

[gradx,grady,gradz] = ComputeGradientMatrix(G);
M=[gradx;grady;gradz];

% compute rotation matrix per face
rot=sparse(3*Nf);
i=1;
Fni=Fn(:,i); % current face normal
rot(i+[0,Nf,2*Nf],i+[0,Nf,2*Nf])=-[0,-Fni(3), Fni(2);Fni(3),0,-Fni(1);-Fni(2),Fni(1),0];
for i=2:Nf
    Fni=Fn(:,i); % current face normal
        %        
    rot(i+[0,Nf,2*Nf],i+[0,Nf,2*Nf])=[0,-Fni(3), Fni(2);Fni(3),0,-Fni(1);-Fni(2),Fni(1),0];
end

A=[M,-rot*M];
%C=sparse(1:2*numel(bIdx),[bIdx;bIdx+Nv],1,2*numel(bIdx),2*Nv); %constrain x and y
C=[sparse(1:numel(bIdx),bIdx,1,numel(bIdx),Nv),zeros(numel(bIdx),Nv)]; % constraint only x
b=[zeros(size(A,2),1);reshape(bValues',[],1)];
x=[A'*A,C';C,zeros(size(C,1))]\b;
x=x(1:size(A,2));
U=[reshape(x,Nv,2)';zeros(1,Nv)];
K=Mesh('VF',U,F);
% K=MeshGraph('VF',U,F);