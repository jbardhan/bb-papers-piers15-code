function [pqrData,a,b,c] = transformPqrToEllipsoidAxes(pqrData)

[a,b,c,X,Y,Z,r0,M] = computeEffectiveEllipsoid(pqrData.xyz,2* ...
					       ones(length(pqrData.R),1));

xyz = pqrData.xyz;

xyz = xyz - ones(length(pqrData.R),1)*r0';
xyz = ([X Y Z]'*xyz')';

pqrData.xyz = xyz;
