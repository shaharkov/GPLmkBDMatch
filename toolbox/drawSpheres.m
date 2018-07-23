function drawSpheres(X,r,c)

if size(c,1)==1
    c = repmat(c,[size(X,1),1]);
end

N = 30;
[x,y,z] = sphere(N);
for ii = 1:size(X,1)
    surf(r*x+X(ii,1), r*y+X(ii,2), r*z+X(ii,3),'edgecolor','none','facecolor',c(ii,:));
end