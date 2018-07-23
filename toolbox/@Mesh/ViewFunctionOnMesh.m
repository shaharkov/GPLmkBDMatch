function h = ViewFunctionOnMesh(G, color_value, options)
%VIEWFUNCTIONONMESH: visualize a function on a triangular mesh
%   This is the full version, the function needs to be defined at all
%   points of the mesh (no interpolation support)
%   The red color stands for positive values, the darker red means the
%   more positive value; the blue color stands for negative value, the
%   darker blue means the more negative value.
%
%   Tingran Gao, Duke University
%   Email: trgao10@math.duke.edu
%   Feb 5, 2015
%

if ~isfield(options, 'mode')
    options.mode = 'rb';
end

if strcmpi(options.mode, 'rb')
    ub_pos = [0 0 1];
    ub_neg = [1 0 0];
    lb = [0.9 0.9 0.8];
    
    pos = zeros(size(color_value));
    pos(color_value>0) = color_value(color_value>0);
    neg = pos-color_value;
    pos = pos/max(pos);
    neg = neg/max(neg);
    
    color_data(color_value>0, :) = pos(color_value>0, :)*ub_pos + (1-pos(color_value>0, :))*lb;
    color_data(color_value<0, :) = neg(color_value<0, :)*ub_neg + (1-neg(color_value<0, :))*lb;
elseif strcmpi(options.mode, 'native')
    color_data = color_value;
end

% figure;
colormap('parula');
h = G.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data, 'CDataMapping', 'scaled', 'EdgeColor', 'none', 'FaceAlpha', 1, 'AmbientStrength',0.3,'SpecularStrength',0.0));
set(gcf, 'ToolBar', 'none');

hold on;
camlight('headlight');
camlight(180,0);
lighting phong;

end

