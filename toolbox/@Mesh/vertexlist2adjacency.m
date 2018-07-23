function A = vertexlist2adjacency(Vlist)
n=numel(Vlist);
A=zeros(n);
for i=1:n
    A(i,Vlist{i})=1;
end

