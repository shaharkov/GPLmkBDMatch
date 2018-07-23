function [WKS,E,PHI,L] = ComputeWKS(G,options)

% function [WKS,E,PHI,L] = compute_WKS(vertices,faces) compute
%   the Wave Kernel Signature of triangle mesh given by
%   [vertices,faces]
%   
%   INPUT:
%   vertices is (number of vertices) x 3 matrix
%   faces is a (number of faces) x 3 matrix
%   
%   OUTPUT:
%   WKS is the (number of vertices) x 100 WKS matrix 
%   E is the vector of LB eigenvalues (by default of size 300 x 1)
%   PHI is the (number of vertices x 300) matrix of LB eigenfunctions 
%   L is the cotan Laplace-Beltrami matrix
%
%   The main parameter to adjust depending on your task is wks_variance


vertices = G.V';
faces = G.F';
n_eigenvalues = getoptions(options,'Nev',300); % number of eigenvalues used for computations
% depending on the application, you can use less than 300
N = getoptions(options,'N',100); % number of evaluations of WKS
wks_variance = getoptions(options,'wks_variance',6); % variance of the WKS gaussian (wih respect to the 
% difference of the two first eigenvalues). For easy or precision tasks 
% (eg. matching with only isometric deformations) you can take it smaller



%% basic quantities
num_vertices = size(vertices,1);
num_faces = size(faces,1);

%% compute cotan LB matrix
% L = G.ComputeLaplacian;

angles = 0*faces;
squared_edge_length = 0*faces;

for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    pp = vertices(faces(:,i2),:) - vertices(faces(:,i1),:);
    qq = vertices(faces(:,i3),:) - vertices(faces(:,i1),:);
    % normalize the vectors
    pp = pp ./ repmat( max(sqrt(sum(pp.^2,2)),eps), [1 3] );
    qq = qq ./ repmat( max(sqrt(sum(qq.^2,2)),eps), [1 3] );
    % compute angles
    angles(:,i1) = acos(sum(pp.*qq,2));
    squared_edge_length(:,i1) = sum((vertices(faces(:,i2)) - vertices(faces(:,i3))).^2,2);
end
clear pp qq;

%then compute L
L = sparse(num_vertices,num_vertices);
for i=1:3
    i1 = mod(i-1,3)+1;
    i2 = mod(i  ,3)+1;
    i3 = mod(i+1,3)+1;
    L = L + sparse(faces(:,i1),faces(:,i2),-cot(angles(:,i3)),...
        num_vertices,num_vertices,num_faces);       
end

L = 1/2 * (L + L');
L = sparse(1:num_vertices,1:num_vertices,-sum(L,2),num_vertices,num_vertices,...
    num_vertices) + L;

fprintf('done. \n');
%% compute Voronoi areas for vertices (following the algorithm decribed in 
%   Discrete Differential-Geometry Operators for Triangulated 2-Manifolds, 
%   Mark Meyer, Mathieu Desbrun, Peter Schroeder and Alan H. Barr, VisMath 2002. 

fprintf('Computing area of vertex Voronoi cells...');

% compute area of triangles
faces_area = zeros(num_faces,1);
for i = 1:3
    faces_area = faces_area +1/4 * (squared_edge_length(:,i).*1./tan(angles(:,i)));
end

% compute area of Voronoi cells
A = zeros(num_vertices,1);
for j = 1:num_vertices
    for i = 1:3
        i1 = mod(i-1,3)+1;
        i2 = mod(i,3)+1;
        i3 = mod(i+1,3)+1;
        ind_j = find(faces(:,i1) == j);
        for l = 1:size(ind_j,1)
            face_index = ind_j(l);
            if (max(angles(face_index,:)) <= pi/2)
                A(j) = A(j) + 1/8 * (1/tan(angles(face_index,i2))* ...
                    squared_edge_length(face_index,i2) + ... 
                    1/tan(angles(face_index,i3))*...
                    squared_edge_length(face_index,i3));
            elseif angles(face_index,i1) > pi/2
                A(j) = A(j) + faces_area(face_index)/2;
            else
                A(j) = A(j) + faces_area(face_index)/4;
            end
        end        
    end
end

A = max(A,1e-8);
area = sum(A);
A = A/area;
Am = sparse([1:length(A)], [1:length(A)], A);

clear angles faces_area;
fprintf('done. \n');

%% compute LB eigenvalues and eigenfunctions 

fprintf('Computing Laplace-Beltrami eigenfunctions...');

% solve generalized eigenvalue problem
options.disp = 0;
[PHI,E] = eigs(L,Am,n_eigenvalues,-1e-5,options);
E = diag(E);
E=abs(real(E));
[E,idx] = sort(E);
PHI = PHI(:,idx);
PHI=real(PHI);

fprintf('done. \n');

%% compute WKS 

fprintf('Computing WKS...');

WKS=zeros(num_vertices,N);

log_E=log(max(abs(E),1e-6))';
e=linspace(log_E(2),(max(log_E))/1.02,N);  
sigma=(e(2)-e(1))*wks_variance;

C = zeros(1,N); %weights used for the normalization of f_E

for i = 1:N
    WKS(:,i) = sum(PHI.^2.*...
        repmat( exp((-(e(i) - log_E).^2) ./ (2*sigma.^2)),num_vertices,1),2);
    C(i) = sum(exp((-(e(i)-log_E).^2)/(2*sigma.^2)));
end

% normalize WKS
WKS(:,:) = WKS(:,:)./repmat(C,num_vertices,1);

fprintf('done. \n');

end

