function L = CORR_compute_laplacian_tension(vertex,face)

%%NEW OPTION - raif
nopts = size(vertex,1);
notrg = size(face,1);
triang = face;
coord = vertex;


% tic;
bir=[1 1 1];
trgarea=zeros(1,notrg);
for i=1:notrg
    trgx=[coord(triang(i,1),1) coord(triang(i,2),1) coord(triang(i,3),1)];
    trgy=[coord(triang(i,1),2) coord(triang(i,2),2) coord(triang(i,3),2)];
    trgz=[coord(triang(i,1),3) coord(triang(i,2),3) coord(triang(i,3),3)];
    aa=[trgx; trgy; bir];
    bb=[trgx; trgz; bir];
    cc=[trgy; trgz; bir];
    area=sqrt(det(aa)^2+det(bb)^2+det(cc)^2)/2;
    trgarea(i)=area;
end


%find the approximate voronoi area of each vertex
AM = zeros(nopts, 1);
for i=1:notrg
    AM(triang(i,1:3)) = AM(triang(i,1:3)) + trgarea(i)/3;
end

% T = sparse([1:nopts], [1:nopts], (AM), nopts, nopts, nopts);
temp = 1./AM;
temp = temp./max(temp);
T = sparse([1:nopts], [1:nopts], temp, nopts, nopts, nopts);




A = sparse(nopts, nopts);

for i=1:notrg
    for ii=1:3
        for jj=(ii+1):3
            kk = 6 - ii - jj; % third vertex no
            v1 = triang(i,ii);
            v2 = triang(i,jj);
            v3 = triang(i,kk);
            e1 = [coord(v1,1) coord(v1,2) coord(v1,3)] - [coord(v2,1) coord(v2,2) coord(v2,3)];
            e2 = [coord(v2,1) coord(v2,2) coord(v2,3)] - [coord(v3,1) coord(v3,2) coord(v3,3)];
            e3 = [coord(v1,1) coord(v1,2) coord(v1,3)] - [coord(v3,1) coord(v3,2) coord(v3,3)];
            cosa = e2* e3'/sqrt(sum(e2.^2)*sum(e3.^2));
            sina = sqrt(1 - cosa^2);
            cota = cosa/sina;
            w = 0.5*cota;
            A(v1, v1) = A(v1, v1) - w;
            A(v1, v2) = A(v1, v2) + w;
            A(v2, v2) = A(v2, v2) - w;
            A(v2, v1) = A(v2, v1) + w;
        end
    end

end
%A = -A;
% L=T*A;
L=A;
return;

%%OLD OPTION - check!
% conformal laplacian

[vertex,face] = CORR_check_face_vertex(vertex,face);

%nface = size(face,1);
n = max(max(face));

L = sparse(n,n);

ring = CORR_compute_vertex_face_ring(face);
bad_row=[];
for i = 1:n
    for b = ring{i}
%         %for isolated vertices
%         if(isempty(b))
%             L(i,i)=1;
%         else
            % b is a face adjacent to a
            bf = face(:,b);
            % compute complementary vertices
            if bf(1)==i
                v = bf(2:3);
            elseif bf(2)==i
                v = bf([1 3]);
            elseif bf(3)==i
                v = bf(1:2);
            else
                error('Problem in face ring.');
            end
            j = v(1); k = v(2);
            vi = vertex(:,i);
            vj = vertex(:,j);
            vk = vertex(:,k);
            % angles
            alpha = myangle(vk-vi,vk-vj);
            beta = myangle(vj-vi,vj-vk);
            % add weight
            cot_alpha = cot(alpha);
            cot_beta = cot(beta);
                L(i,j) = L(i,j) +cot_alpha;
                L(i,k) = L(i,k) + cot_beta;
            
%         end
    end
end


L = L - diag(sum(L,2));


return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );