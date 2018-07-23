function [] = CORR_draw_spheres(P,R,color)
% h- axes handles
% P nx3 locations
% R nx1 radius
% color nx3 colors

%if only on color make it for every point
if(length(R) == 1)
    R = R*ones(size(P,1),1);
end

if(size(color,1) == 1)
    color = ones(size(P,1),1)*color;
end

[x, y, z] = sphere(30);

for k=1:size(P,1)
    xx = R(k)*x+ P(k,1);
    yy = R(k)*y+ P(k,2);
    zz = R(k)*z+ P(k,3);
    mesh(xx,yy,zz, 'EdgeColor','none','FaceColor', color(k,:),'AmbientStrength',0.5);
end

% lighting phong
% material shiny;

