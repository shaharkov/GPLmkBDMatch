function [Gx,Gy,Gz] = ComputeGradientMatrix(G)

%%% compute gradient of domainMesh on each face
XF = @(i)G.V(:,G.F(i,:));
Na = cross( XF(2)-XF(1), XF(1)-XF(3) );
amplitude = @(X)sqrt( sum( X.^2 ) );
A = amplitude(Na)/2;
normalize = @(X)X ./ repmat(amplitude(X), [3 1]);
N = normalize(Na);
I = []; J = []; V = []; % indexes to build the sparse matrices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % vector N_f^e_i
    wi = cross(XF(t)-XF(s),N);
    % update the index listing
    I = [I, 1:G.nF];
    J = [J, G.F(i,:)];
    V = [V, wi];
end

dA = spdiags(1./(2*A(:)), 0, G.nF, G.nF);

Gx = dA*sparse(I,J,V(1,:),G.nF,G.nV);
Gy = dA*sparse(I,J,V(2,:),G.nF,G.nV);
Gz = dA*sparse(I,J,V(3,:),G.nF,G.nV);

end

