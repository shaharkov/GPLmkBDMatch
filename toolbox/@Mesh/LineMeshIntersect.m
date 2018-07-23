function [pt] = LineMeshIntersect(G, segment)
%LINEMESHINTERSECT: return the intersection of a line segment with the mesh
%   
%   'segment' is a 2x3 matrix, representing the beginning/end point
%
%   Tingran Gao, Duke University
%   Email: trgao10@math.duke.edu
%   Feb 6, 2015
%

pt = [];
frontIntersect = 0;

for j=1:G.nF
    TriV = G.V(:,G.F(:,j));
    t = [(segment(2,:)-segment(1,:))',TriV(:,2)-TriV(:,1),TriV(:,3)-TriV(:,1)]\(segment(2,:)'-TriV(:,1));
    if ((t(1) < 0 || t(1) > 1) || (t(2) < 0 || t(2) > 1) || (t(3) < 0 || t(3) > 1))
        continue;
    end
    if (t(1) > frontIntersect)
        frontIntersect = t(1);
        pt = segment'*[t(1);1-t(1)];
    end
end

% hold on
% scatter3(segment(1,1),segment(1,2),segment(1,3),'r','filled');
% scatter3(segment(2,1),segment(2,2),segment(2,3),'b','filled');
% scatter3(pt(1), pt(2), pt(3), 20, 'g', 'filled');
% hold off


end

