function K = ComputeCPMS(G)
V=G.V';
F=G.F';
nv = size(V,1);
nf = size(F,1);

% Find orig edge lengths and angles
EL = G.ComputeEdgeLengths;
A = G.ComputeTriangleAngles(EL);

% Build adjacency and face matrices
I = [F(:,1);F(:,2);F(:,3)];
J = [F(:,2);F(:,3);F(:,1)];
S = [1:nf,1:nf,1:nf];
E = sparse(I,J,S',nv,nv);
D = double(E > 0);
D = D + D';
D = double(D > 0);

% Find boundary vertices
[I,J] = find(E);
BV = sparse(nv,1);
for m=1:length(I)
    i = I(m); j = J(m);
    if xor(E(i,j),E(j,i))
        BV(i) = 1;
        BV(j) = 1;
    end
end


% Find original curvature - Korig
I = reshape(F,nf*3,1);
J = ones(nf*3,1);
S = reshape(A,nf*3,1);
SA = sparse(I,J,S,nv,1);    
SA = full(SA);
Korig = 2*pi*ones(nv,1) - pi*BV - SA;


% Find Laplacian - L
I = [F(:,1);F(:,2);F(:,3)];
J = [F(:,2);F(:,3);F(:,1)];
S = 0.5*cot([A(:,3);A(:,1);A(:,2)]);
In = [I;J;I;J];
Jn = [J;I;I;J];
Sn = [-S;-S;S;S];
L = sparse(In,Jn,Sn,nv,nv);


% Find PHI (1)
locs = find(BV==0);
PHI = zeros(nv,1);
PHI(locs) = -L(locs,locs)\Korig(locs);


% Find new metric - ELnew (2)
S1 = scale_phi(PHI,F(:,2),F(:,3));
S2 = scale_phi(PHI,F(:,3),F(:,1));
S3 = scale_phi(PHI,F(:,1),F(:,2));
S = [S1, S2, S3];
ELnew = EL.*S;



% Embed the metric
V = full(lscm(ELnew,F));
K=Mesh('VF',V',F');
end

function S = scale_phi(PHI, I, J)
S = exp(PHI(J))+exp(PHI(I));
end

function V = lscm(ELf,T)

nv = max(max(T));
nf = length(T);

L1 = ELf(:,1)'; L2 = ELf(:,2)'; L3 = ELf(:,3)';

A1 = real(acos((L3.^2+L2.^2-L1.^2)./(2.*L3.*L2)));
WXX = L2./L3.*cos(A1);
WXY = L2./L3.*(-sin(A1));
neq = 2*nf;
I = reshape(repmat([1:neq]',1,6)',6*neq,1);
Ti = reshape(repmat([1:nf]',1,4)',4*nf,1);
TT = T(Ti,:);
TT(2:2:end) = TT(2:2:end)+nv;
J = reshape(TT',1,2*neq*3)';
S = [WXX-1; -WXX; ones(1,nf); WXY; -WXY; zeros(1,nf); ...
     -WXY;   WXY; zeros(1,nf); WXX-1; -WXX; ones(1,nf)];
S = reshape(S,size(S,1)*size(S,2),1);
A = sparse(I,J,S,2*nf+4,2*nv); 


v1 = T(nf,1); v2 = T(nf,2);

A(end-3,v1) = 1;
A(end-2,v2) = 1;
A(end-1,nv+v1) = 1;
A(end,nv+v2) = 1;
B = sparse(2*nf+4,1);
B(end-1) = ELf(nf,3);

X2Dl = A\B;

V = [X2Dl(1:nv,1),X2Dl(nv+1:end,1)];

s = mean(max(V)-min(V));
V = V./s*10;
end
