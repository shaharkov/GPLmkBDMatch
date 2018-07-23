function [K,V2V,E2E] = Soupalize(G)
%SOUPALIZE returns a polygon soup
%   K: Returned soup
%   V2V: An "adjacency" #Vx3#F matrix between vertices of G and K
%   E2E: An "adjacency" #Ex3#F matrix between edges of G and K
if iscell(G.F)
    error('Not implemented for non-triangular meshes yet');
end
if ~isempty(G.Soup)
    K=G.Soup;
    V2V=K.V2V;
    E2E=K.E2E;
    return;
end
V=G.V; Nv=size(V,2);
F=G.F; Nf=size(F,2);


U=zeros(3,3*Nf);
U(:,1:3*Nf)=V(:,F(:));
K=Mesh('VF',U,reshape(1:3*Nf,3,[]));
if nargout>1 
    V2V=sparse(F(:),1:3*Nf,1,Nv,3*Nf);
    K.V2V=V2V;
end
if nargout>2
    Ic_G = G.ComputeV2E;
    Ic_K = K.ComputeV2E;
    % This find the parent j of edge i in the soup
    E2E=mod(find(abs(Ic_G'*V2V*Ic_K)==2),size(Ic_G,2)); E2E(E2E==0)=size(Ic_G,2);
    E2E=sparse(1:size(Ic_K,2),E2E,1,size(Ic_K,2),size(Ic_G,2)); E2E=E2E';
    % index of last non-zero element of each row of E2E
    %j = sum(cumsum((E2E(:,end:-1:1) ~= 0), 2) ~= 0, 2);
    j=zeros(size(E2E,1),1);
    for i=1:size(E2E,1)
        j(i)= find(E2E(i,:), 1, 'last');
    end
    E2E(sub2ind(size(E2E),1:size(E2E,1),j'))=-1;
    K.E2E=E2E;
end

G.Soup=K;

end

