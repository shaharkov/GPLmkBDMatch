function Write(G,filename,format,options)

if nargin < 4
    options = struct('pointCloud', 0);
end
% options.pointCloud = getoptions(options, 'pointCloud', 0);

switch format
    case 'off'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end
        
        % header
        fprintf(fid, 'OFF\n');
        if options.pointCloud==0
            fprintf(fid, '%d %d 0\n', length(G.V), length(G.F));
        else
            fprintf(fid, '%d 0 0\n', length(G.V));
        end
        
        % write the points & faces
        fprintf(fid, '%f %f %f\n', G.V);
        if options.pointCloud==0
            fprintf(fid, '3 %d %d %d\n', G.F-1);
        end
        
        fclose(fid);
    case 'obj'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end
        
        % vertex coordinates
        fprintf(fid, 'v %f %f %f\n', G.V);
        
        % Texture coordinates
        if isfield(options, 'Texture')
            fprintf(fid, 'vt %f %f\n', options.Texture.Coordinates(1:2,:));
            fprintf(fid, 'f %d/%d %d/%d %d/%d\n', kron(G.F',[1,1])');
        else
            fprintf(fid, 'f %d %d %d\n', G.F');
        end
        fclose(fid);
end