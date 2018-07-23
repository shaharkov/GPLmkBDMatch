function TEETH_output_vtk_texture(filename,V1,F1,UV)

ofid = fopen(filename,'w');
fprintf(ofid, '# vtk DataFile Version 3.0\n');
fprintf(ofid,'vtk output\n');
fprintf(ofid,'ASCII\n');
fprintf(ofid,'DATASET POLYDATA\n');
fprintf(ofid,'POINTS %d float\n', size(V1,1));
fprintf(ofid,'%g %g %g\n', V1');
fprintf(ofid,'POLYGONS %d %d\n', size(F1,1), 4*size(F1,1));
fprintf(ofid,'3 %d %d %d\n', F1'-1);
fprintf(ofid,'\n');

%texture coordinates
UV(:,1) = 0.5*(UV(:,1)+1);
UV(:,2) = 0.5*(UV(:,2)+1);

fprintf(ofid,'POINT_DATA %d\n', size(V1,1));
fprintf(ofid,'TEXTURE_COORDINATES TCoords 2 float\n');
fprintf(ofid,'%g %g\n', UV');

fclose(ofid);
