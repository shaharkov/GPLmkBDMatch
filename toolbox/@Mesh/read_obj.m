  function [V,F,Fs] = read_obj( filename )
% function [V,F] = read_obj( filename )
%
% Loads a mesh in Wavefront OBJ format.

   fp = fopen( filename, 'r' );
   if( fp == -1 )
      disp( sprintf( 'Error: could not read mesh file "%s"\n', filename ));
      return;
   end

   V = zeros( 3, 0 );
   F = zeros( 3, 0 );
   while( ~feof( fp ))
       line = fscanf( fp, '%c %f %f %f' );
       label = line( 1 );
       values = line( 2:4 );

       if( label == 'v' ) % vertex
          V = [ V, values(1:3) ];
       end

       if( label == 'f' ) % face
          F = [ F, values ];
       end
   end
   fclose( fp );
   Fs=repmat(3,1,size(F,2));
end

