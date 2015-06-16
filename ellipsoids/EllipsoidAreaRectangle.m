function [ AreaEllipsoid ] = EllipsoidAreaRectangle(xc,yc,zc,xr,yr,zr,n)
% EllipsoidAreaRectangle: outputs the Area of an Ellipsoid using n
%   rectangles per side.
%
%   AreaEllipsoid = EllipsoidAreaRectangle(xc, yc, zc, xr, yr, zr, n)
%       Inputs: (xc,yc,zc) - center of ellipsoid
%               (xz,yr,zr) - x, y, z semi-axes
%               n          - number of rectangles in longitude and latitude
%       Output: AreaEllipsoid
%
% Written by Noah Wieckowski and Danish Tejani, 2014.

[x,y,z] = ellipsoid(xc,yc,zc,xr,yr,zr,n);

RectCounter = 1;

for zIndex = 1:n
    
    for xIndex = 1:n
        
        v1 = [x(zIndex, xIndex) y(zIndex,xIndex) z(zIndex, xIndex)];
        v2 = [x(zIndex, xIndex+1) y(zIndex, xIndex+1) z(zIndex, xIndex+1)];
        v3 = [x(zIndex+1, xIndex+1) y(zIndex+1, xIndex+1) z(zIndex+1, xIndex+1)];
        v4 = [x(zIndex+1, xIndex) y(zIndex+1, xIndex) z(zIndex+1, xIndex)];
        
        Rect(RectCounter,:) = [v1 v2 v3 v4];
        
        Area(RectCounter) = 0.5*norm(cross((v4-v1),(v2-v1))) + 0.5*norm(cross((v2-v3),(v4-v3)));
        
        RectCounter = RectCounter+1;
    end
    
end

AreaEllipsoid = sum(Area);

end