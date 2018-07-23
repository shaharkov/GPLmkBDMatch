function [M, edge_to_vertex_map, edge_num] = CORR_compute_edge_numbering(F)

num_of_vertices = max(max(F));

%M=-ones(num_of_vertices,num_of_vertices);
M = sparse(num_of_vertices, num_of_vertices);

cnt=0;
for i=1:size(F,1)
    n1 = F(i,1);
    n2 = F(i,2);
    n3 = F(i,3);
    
    if(  M(n1,n2) == 0 )
        cnt=cnt+1;
        M(n1,n2)=cnt;
        M(n2,n1)=cnt;
    end
    
    if( M(n3,n1) == 0)
        cnt=cnt+1;
        M(n1,n3)=cnt;
        M(n3,n1)=cnt;        
    end
    
    if( M(n2,n3) == 0)
        cnt=cnt+1;
        M(n3,n2)=cnt;
        M(n2,n3)=cnt;        
    end
   
    
    
end

edge_num = cnt;

edge_to_vertex_map = zeros(edge_num,2);
ind = find(M ~= 0);
dimM=size(M,1);
for k=1:length(ind)
    %where ind is in the matrix M\
    temp = mod(ind(k),dimM);
    edge_to_vertex_map(M(ind(k)),:) = [int32(1 + ( (ind(k)-temp)/dimM ) ) , temp ];
    
end