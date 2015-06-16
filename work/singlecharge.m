visualizeOn = 0;
addpath('../pointbem');
addpath('../panelbem');
addpath('../testasymmetry');
addpath('../ellipsoids');


standardAtomSize = 1.5;
nonlinearResponseDistance = 5;
deltaRadius = 0.5;

%%%%%%%%%  SOLVENT MODEL
loadConstants;
origin = [0 0 0];
epsIn = 1;
epsOut = 80;
kappa = 0.0; % we're not treating salt in this work. change to
             % 0.125 to model physiological salt conditions
conv_factor = 332.112; % converts to kcal/mol
staticpotential = 10.7; % obtained in earlier work
sternWidth = 2.0; % Angstroms
symParams  = struct('alpha',0.0, 'beta',   0.0,'EfieldOffset', 0.0);
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);



%%%%%%%%% SOLUTE MODEL
density = 1;

minRadius   = 3;
deltaRadius = 1;
maxRadius   = 10;
atomRadius = 1.5;

minQ = -1;
maxQ = +1;
numQ = 2;

radius = minRadius:deltaRadius:maxRadius;
qList = minQ:(maxQ-minQ)/(numQ-1):maxQ;

for i=1:length(radius)
  curRadius = radius(i);

  numDielPoints(i)  = ceil(4 * pi * density * curRadius^2);
  dielSurfData   = makeSphereSurface(origin, curRadius, numDielPoints(i));
  numSternPoints(i) = ceil(4 * pi * density * (curRadius+sternWidth)^2);
  sternSurfData  = makeSphereSurface(origin, curRadius+sternWidth, numSternPoints(i));
  numTotalPoints(i) = numDielPoints(i) + numSternPoints(i);

  for j=1:length(qList)
    curQ = qList(j);
    pqrBuried  = struct('xyz',[0 0 0],'q',curQ,'R',0);
    pqrSurface = struct('xyz',[curRadius-atomRadius 0 0],'q',curQ,'R',0);

    bemEcfAsym      = makeBemEcfQualMatrices(dielSurfData, pqrBuried,  epsIn, epsOut);
    bemStern        = makeBemSternMatrices(dielSurfData, sternSurfData, ...
					   pqrBuried, epsIn, epsOut, kappa);
    [phiReacPointSternAsym, phiBndy,dPhiBndy] = ...
	solveConsistentSternAsym(dielSurfData, sternSurfData, pqrBuried, ...
				 bemStern, epsIn, epsOut, kappa, ...
				 conv_factor, asymParams); % symParams!!
    [phiReacPointSternSym, phiBndy,dPhiBndy] = ...
	solveConsistentSternAsym(dielSurfData, sternSurfData, pqrBuried, ...
				 bemStern, epsIn, epsOut, kappa, ...
				 conv_factor, symParams); % symParams!!
    dG_asym_buried(i,j) = 0.5 * pqrBuried.q'*phiReacPointSternAsym;
    dG_sym_buried(i,j) = 0.5 * pqrBuried.q'*phiReacPointSternSym;
    withstatic_asym_buried(i,j) = dG_asym_buried(i,j)+sum(pqrBuried.q)*staticpotential;
    withstatic_sym_buried(i,j) = dG_sym_buried(i,j)+sum(pqrBuried.q)*staticpotential;
    
    bemEcfAsym      = makeBemEcfQualMatrices(dielSurfData, pqrSurface,  epsIn, epsOut);
    bemStern        = makeBemSternMatrices(dielSurfData, sternSurfData, ...
					   pqrSurface, epsIn, epsOut, kappa);
    [phiReacPointSternAsym, phiBndy,dPhiBndy] = ...
	solveConsistentSternAsym(dielSurfData, sternSurfData, pqrSurface, ...
				 bemStern, epsIn, epsOut, kappa, ...
				 conv_factor, asymParams); % symParams!!
    [phiReacPointSternSym, phiBndy,dPhiBndy] = ...
	solveConsistentSternAsym(dielSurfData, sternSurfData, pqrSurface, ...
				 bemStern, epsIn, epsOut, kappa, ...
				 conv_factor, symParams); % symParams!!
    dG_asym_surface(i,j) = 0.5 * pqrSurface.q'*phiReacPointSternAsym;
    dG_sym_surface(i,j) = 0.5 * pqrSurface.q'*phiReacPointSternSym;
    withstatic_asym_surface(i,j) = dG_asym_surface(i,j)+sum(pqrSurface.q)*staticpotential;
    withstatic_sym_surface(i,j) = dG_sym_surface(i,j)+sum(pqrSurface.q)*staticpotential;
  end
end