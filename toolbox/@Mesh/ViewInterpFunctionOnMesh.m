function ViewInterpFunctionOnMesh(GM, color_value, options)
%VIEW_MESH_FUNCTION_ON_SAMPLES Summary of this function goes here
%   Detailed explanation goes here

if (~isempty(findobj('Type','figure')))
    camUpVector = get(gca, 'CameraUpVector');
    camPosition = get(gca, 'CameraPosition');
    camTarget = get(gca, 'CameraTarget');
    camViewAngle = get(gca, 'CameraViewAngle');
end

figure;
set(gcf, 'ToolBar', 'none');

interp_color_ub_pos = [0 0 1];
interp_color_ub_neg = [1 0 0];
interp_color_lb = [0.9 0.9 0.8];

color_data = repmat(interp_color_lb', 1, max(size(GM.V)));
invisible = ones(max(size(GM.V)),1);
SampNumb = GM.Aux.SampNumb;
SampInd = GM.Aux.VertSampInd;
color_value_pos = zeros(size(color_value));
color_value_pos(color_value>0) = color_value(color_value>0);
color_value_neg = color_value_pos-color_value;
color_value_pos = color_value_pos/max(color_value_pos);
color_value_neg = color_value_neg/max(color_value_neg);

if strcmpi(options.boundary,'off')
    bvs_inds = compute_boundary(GM.F);
    [D,~,~] = PerformFastMarching(GM,bvs_inds);
    SampToBvs = D(SampInd);
    SampToBvs(isinf(SampToBvs)) = 0;
    
    color_value_pos = color_value_pos.*SampToBvs;
    color_value_pos = color_value_pos/max(color_value_pos);
    
    color_value_neg = color_value_neg.*SampToBvs;
    color_value_neg = color_value_neg/max(color_value_neg);
end

% interpolate positive values
for j=1:SampNumb
    color_float_index = color_value_pos(j);
    if (color_float_index>0)
        color_position = find(GM.Aux.V2S == SampInd(j));
        % larger color_float_index gets more interp_ub_color
        interp_color = interp_color_ub_pos'*color_float_index + interp_color_lb'*(1-color_float_index);
        color_data(:, color_position) = repmat(interp_color, 1, length(color_position));
        
        if strcmpi(options.transparency, 'positive')
            invisible(color_position) = 0;
        end
    end
end

% interpolate negative values
for j=1:SampNumb
    color_float_index = color_value_neg(j);
    if (color_float_index>0)
        color_position = find(GM.Aux.V2S == SampInd(j));
        % larger color_float_index gets more interp_ub_color
        interp_color = interp_color_ub_neg'*color_float_index + interp_color_lb'*(1-color_float_index);
        color_data(:, color_position) = repmat(interp_color, 1, length(color_position));
        
        if strcmpi(options.transparency, 'negative')
            invisible(color_position) = 0;
        end
    end
end

if strcmpi(options.transparency, 'none')
    GM.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data', 'CDataMapping', 'scaled', 'EdgeColor', 'none', 'FaceVertexAlphaData', invisible, 'FaceAlpha', 1, 'AmbientStrength',0.3,'SpecularStrength',0.0));
else
    GM.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data', 'CDataMapping', 'scaled', 'EdgeColor', 'none', 'FaceVertexAlphaData', invisible, 'FaceAlpha', 'interp', 'AmbientStrength',0.3,'SpecularStrength',0.0));
end
% GM.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data', 'CDataMapping','scaled', 'EdgeColor', 'none', 'FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
hold on;
camlight('headlight');
camlight(180,0);
lighting phong;
title(GM.Aux.name);

if (exist('camUpVector', 'var'))
    set(gca, 'CameraUpVector', camUpVector);
    set(gca, 'CameraPosition', camPosition);
    set(gca, 'CameraTarget', camTarget);
    set(gca, 'CameraViewAngle', camViewAngle);
end

end

