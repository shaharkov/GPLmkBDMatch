function Fn=ComputeFaceNormals(G)
F=G.F;
V=G.V;
nf=size(F,2);
Fn=zeros(3,nf);
if size(F,1)~=3
    error('Not a triangular mesh!');
end
for i=1:nf
    e1=V(:,F(2,i))-V(:,F(1,i));
    e2=V(:,F(3,i))-V(:,F(1,i));
    Fn(:,i)=cross(e1,e2);
    Fn(:,i)=Fn(:,i)./sqrt(sum(Fn(:,i).^2));
end
