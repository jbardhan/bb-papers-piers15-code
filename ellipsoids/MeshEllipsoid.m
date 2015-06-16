function [surfData,AreaEllipsoid,TriAreaEllipsoid] = MeshEllipsoid(a, b, c, n)
% MeshEllipsoid: returns a triangulated ellipsoid using n
%   rectangles per side.
%
%   [surfData,AreaEllipsoid,TriArea] = MeshEllipsoid(a, b, c, n)
%
%       Inputs: (a, b, c) - semi-axes (will be along x, y, z axes, respectively)
%               n          - number of rectangles in longitude and latitude
%
%       Output: surfData         - a data structure compatible with our panel BEM
%               AreaEllipsoid    - the total area of the mesh,
%                      computed with rectangles
%               TriAreaEllipsoid - the total area of the mesh,
%                      computed with rectangles
%
% Based on EllipsoidAreaRectangle (Noah Wieckowski and Danish
% Tejani, 2014). Modifications begun by Jay Bardhan, 2015.

xc = 0; yc = 0; zc = 0;
[x,y,z] = ellipsoid(xc, yc, zc, a, b, c, n);

RectCounter = 1;
TriCounter  = 1;

for zIndex = 1:n
    
    for xIndex = 1:n
        
        v1 = [x(zIndex, xIndex) y(zIndex,xIndex) z(zIndex, xIndex)];
        v2 = [x(zIndex, xIndex+1) y(zIndex, xIndex+1) z(zIndex, xIndex+1)];
        v3 = [x(zIndex+1, xIndex+1) y(zIndex+1, xIndex+1) z(zIndex+1, xIndex+1)];
        v4 = [x(zIndex+1, xIndex) y(zIndex+1, xIndex) z(zIndex+1, xIndex)];
        
        Rect(RectCounter,:) = [v1 v2 v3 v4];
        RectArea(RectCounter) = 0.5*norm(cross((v4-v1),(v2-v1))) + 0.5*norm(cross((v2-v3),(v4-v3)));
        RectCounter = RectCounter+1;
        
	% the top and bottom z sections are triangles
	if zIndex > 1
	  Tri(TriCounter,:) = [v1 v2 v4];
	  TriCentroids(TriCounter,:) = mean([v1; v2; v4]);
	  TriNormals(TriCounter,:) = -1 * cross((v4-v1),(v2-v1)); % so normal points out.
	  % note that this is the same as (cross(v2-v1),(v4-v1))!
	  TriArea(TriCounter) = 0.5*norm(TriNormals(TriCounter,:));
	  TriNormals(TriCounter,:) = TriNormals(TriCounter,:)/norm(TriNormals(TriCounter,:));
	  TriCounter = TriCounter + 1;
	end
	
	if zIndex < n
	  Tri(TriCounter,:) = [v2 v3 v4];
	  TriCentroids(TriCounter,:) = mean([v2; v3; v4]);
	  TriNormals(TriCounter,:) = -1 * cross((v2-v3),(v4-v3));
	  TriArea(TriCounter) =  0.5*norm(TriNormals(TriCounter,:));
	  TriNormals(TriCounter,:) = TriNormals(TriCounter,:)/norm(TriNormals(TriCounter,:));
	  TriCounter = TriCounter + 1;
	end

    end
    
end
mesh = Tri;
X = mesh(:, [1 4 7]);
Y = mesh(:, [2 5 8]);
Z = mesh(:, [3 6 9]);
meshData = struct('vert',0,'face',ones(length(TriArea),1), 'X',X','Y',Y','Z',Z', 'normals', TriNormals');

AreaEllipsoid = sum(RectArea);
TriAreaEllipsoid = sum(TriArea);

surfData = struct('meshData',meshData,'areas',TriArea,'centroids',TriCentroids,...
		  'normals',TriNormals);