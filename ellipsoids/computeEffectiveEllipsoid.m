function [a, b, c, X, Y, Z, r0, M] = ...
    computeEffectiveEllipsoid(xyz, R)

% computeEffectiveEllipsoid uses the method of Sigalov, Fenley, and
%   Onufriev (J. Chem. Phys. v.124:124902 (2006) ) to compute the
%   primary axes and semi-axes for a protein whose shape is
%   specified as a set of spheres
%
% [a,b,c,X,Y,Z,r0,M] = computeEffectiveEllipsoid(xyz, R)
%   Inputs:  xyz = [x1 y1 z1; x2 y2 z2; ...];
%              R = vector of atom radii
%
%   Outputs: a, b, c the semi-axis lengths
%            X, Y, Z the semi-axes themselves, column vectors
%            r0 the center of "mass"
%            M  the total "mass" 
%    to go from "global coordinates" (xyzr(1,:) for example)
%       use localXYZ = [X Y Z]' * (xyzr(1,:)-r0)'

numAtoms = size(xyz,1); % size(A,1) gives number of rows
xyz = xyz'; % now it's a set of column vectors. much easier to use! 

% step 1: compute total "mass" = sum of volumes
m = R.^3; 
M = sum(m);

% step 2: compute center of "mass"
r0 = zeros(3,1);
for i=1:numAtoms
  r0 = r0 + m(i) * xyz(:,i);
end
r0 = r0 / M;

% step 3: translate molecule so the center of "mass" is at the
% origin
xyz = xyz - r0 * ones(1,numAtoms);

% step 4: build second moment matrix Imoment and find eigendecomposition
Imoment = zeros(3); % this creates a square matrix of size 3 with all
                    % entries zero
x = xyz(1,:);
y = xyz(2,:);
z = xyz(3,:);
for i=1:numAtoms
  Imoment(1,1) = Imoment(1,1)+m(i)*(y(i)^2+z(i)^2+(2/5)*R(i)^2);
  Imoment(2,2) = Imoment(2,2)+m(i)*(x(i)^2+z(i)^2+(2/5)*R(i)^2);
  Imoment(3,3) = Imoment(3,3)+m(i)*(x(i)^2+y(i)^2+(2/5)*R(i)^2);

  Imoment(1,2) = Imoment(1,2)-m(i)*x(i)*y(i);
  Imoment(1,3) = Imoment(1,3)-m(i)*x(i)*z(i);
  Imoment(2,3) = Imoment(2,3)-m(i)*y(i)*z(i);
end
Imoment(2,1) = Imoment(1,2);
Imoment(3,1) = Imoment(1,3);
Imoment(3,2) = Imoment(2,3);

[V,D] = eig(Imoment);

% step N-1: assign Ixx,Iyy,Izz to increasing eigenvalues, and
% associate each with the appropriate eigenvector;
[lambda,I] = sort(diag(D),'ascend');
Ixx = lambda(1); X = V(:,I(1));
Iyy = lambda(2); Y = V(:,I(2));
Izz = lambda(3); Z = V(:,I(3));

% step N: Eq 17 in Sigalov06
a = sqrt((5/(2*M))*(-Ixx+Iyy+Izz));
b = sqrt((5/(2*M))*(+Ixx-Iyy+Izz));
c = sqrt((5/(2*M))*(+Ixx+Iyy-Izz));