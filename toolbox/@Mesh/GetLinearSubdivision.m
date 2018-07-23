function K =  GetLinearSubdivision(G)
% Linear subdivision for triangle meshes
%
%  Dimensions:
%    vertices: 3xn Vertices
%    faces:    3xn Faces
%  
%  Author: Jesus Mena
%
% maintain ancestor list (which face each face originated from.
% use input ancestor list to decide ancestor of new faces.

	global edgeVertex;
    global newIndexOfVertices;
	vertices=G.V;
    faces=G.F;
    ancestor=G.ancestors;
    newFaces = [];
	newVertices = vertices;
	nVertices = size(vertices,2);
	nFaces    = size(faces,2);
	edgeVertex= zeros(nVertices, nVertices);
	newIndexOfVertices = nVertices;

    %old faces keep their old ancestors
    if isempty(ancestor)
        ancestor=1:nFaces;
    end
    new_ancestor = zeros(4*length(ancestor),1);
    
    % ------------------------------------------------------------------------ %
	% create a matrix of edge-vertices and a new triangulation (newFaces).
    % 
    % * edgeVertex(x,y): index of the new vertex between (x,y)
    %
    %  0riginal vertices: va, vb, vc.
    %  New vertices: vp, vq, vr.
    %
    %      vb                   vb             
    %     / \                  /  \ 
    %    /   \                vp--vq
    %   /     \              / \  / \
    % va ----- vc   ->     va-- vr --vc 
	%
    cntyl=-3;
    for i=1:nFaces
        [vaIndex, vbIndex, vcIndex] = deal(faces(1,i), faces(2,i), faces(3,i));
        
        vpIndex = addEdgeVertex(vaIndex, vbIndex);
        vqIndex = addEdgeVertex(vbIndex, vcIndex);
        vrIndex = addEdgeVertex(vaIndex, vcIndex);
        
        fourFaces = [vaIndex,vpIndex,vrIndex; vpIndex,vbIndex,vqIndex; vrIndex,vqIndex,vcIndex; vrIndex,vpIndex,vqIndex]';
        newFaces  = [newFaces, fourFaces];
        
        %update ancestor list: new 4 faces point to the ancestor of this
        %i-th face
        cntyl=cntyl+4;
        new_ancestor(cntyl:(cntyl+3)) = ancestor(i)*ones(4,1);
        %         new_ancestor = [new_ancestor; ancestor(i)*ones(4,1) ];
    end
    	
    % ------------------------------------------------------------------------ %
	% positions of the new vertices
	for v1=1:nVertices-1
		for v2=v1:nVertices
			vNIndex = edgeVertex(v1,v2);
            if (vNIndex~=0)
 				newVertices(:,vNIndex) = 1/2*(vertices(:,v1)+vertices(:,v2));
            end
        end
    end
 	K=MeshGraph('VF',newVertices,newFaces);
    K.ancestors=new_ancestor;
end

% ---------------------------------------------------------------------------- %
function vNIndex = addEdgeVertex(v1Index, v2Index)
	global edgeVertex;
	global newIndexOfVertices;

	if (v1Index>v2Index) % setting: v1 <= v2
		vTmp = v1Index;
		v1Index = v2Index;
		v2Index = vTmp;
	end;
	
	if (edgeVertex(v1Index, v2Index)==0)  % new vertex
		newIndexOfVertices = newIndexOfVertices+1;
		edgeVertex(v1Index, v2Index) = newIndexOfVertices;
	end;

	vNIndex = edgeVertex(v1Index, v2Index);

    return;
end
