function n=ComputeNonPlanarFaceNormals(G)
nf=size(G.FVIc,2);
n=zeros(3,nf);
%d=zeros(nf,1);
%s=zeros(nf,1);

for k=1:nf
    Q=G.V(:,find(G.FVIc(:,k)))';
    M=sum(Q,1)/size(Q,1);
    Qc=Q-M(ones(size(Q,1),1),:);
    [~,~,v]=svd(Qc'*Qc);
    n(:,k)=v(:,3)';
%    s(k)=ss(3,3);
%    d(k)=-M*v(:,3);
end

end