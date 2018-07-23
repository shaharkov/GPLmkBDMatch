function [OBV,OBE]  = FindOrientedBoundaries( obj )
% FINDORIENTEDBOUNDARIESMATLAB returns a cell array of oriented boundary 
%  vertices (a cell for every boundary) and a cell array for boundary edges
%  a cell for every boundary.
%   input: TriangleMesh object
%   output: OBV cell array with a cell for every boundary vertices
%   output: OBE cell array with a cell for every boundary edges
%
% Created by Nave Zelzer on may 23 2014.

[BV,BE] = obj.FindBoundaries();
if isempty(BE)
    OBV = {};
    OBE = {};
    return;
end
OBV = cell(1,size(BE,2));
OBE = cell(1,size(BE,2));

edge = BE(:,1);
BE(:,1) = [];
to = edge(1);
nb = 1; % number of boundaries
while ~isempty(BE)
    OBV{nb} = [to, OBV{nb}];
    OBE{nb} = [edge, OBE{nb}];
	next = BE(:,(BE(1,:) == to | BE(2,:) == to));
	next(:,ismember(next',edge','rows')) = [];
    if isempty(next)
        edge = BE(:,1);
        to = edge(1);
        BE(:,1) = [];
        nb = nb+1;
        continue;
    else
        next = next(:,1);
    end
	edge = next;
    to = edge(edge ~= to);
	BE(:,ismember(BE',edge','rows')) = [];
end
OBV{nb} = [to, OBV{nb}];
OBE{nb} = [edge, OBE{nb}];
OBV = OBV(~cellfun('isempty',OBE));
OBE = OBE(~cellfun('isempty',OBE));
OBV = OBV{:};
OBE = OBE{:};
end

