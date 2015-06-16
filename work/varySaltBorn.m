addpath('../pointbem');
addpath('../panelbem');
addpath('../testasymmetry');
loadConstants

origin = [0 0 0];
R_list = linspace(1,2.3,20);
% sodium RminOver2 1.41075 is from the newest Roux paper
% chloride RminOver2 2.27 

sternLayerThickness = 2.0;
chargeLocation = [0 0 0];
q_list = [-1 1];
kappa_list = [0.000];

epsIn  =  1;
epsOut = 80;
conv_factor = 332.112;
asymParams = struct('alpha',0.5, 'beta', -60.0,'EfieldOffset',-0.5);
symParams  = struct('alpha',0.0, 'beta', -60.0,'EfieldOffset',-0.5);

k = 1; 
densityDiel = 4;
densityStern = 1;
staticPotential = 10.7; % approximate static potential. weakly size
                        % dependent, see Ashbaugh00, Bardhan12
			
for i=1:length(R_list)
  R = R_list(i);
  Rstern = R + sternLayerThickness;

  numDielPoints(i) = ceil(4 * pi *  densityDiel * R^2);
  numSternPoints(i) = ceil(4 * pi * densityStern * Rstern^2);
  numTotalPoints(i) = numDielPoints(i) + numSternPoints(i);
  
  dielSurfData = makeSphereSurface(origin, R, numDielPoints(i));
  sternSurfData = makeSphereSurface(origin, Rstern, numSternPoints(i));

  for j=1:length(q_list)
    q = q_list(j);
    pqrData = struct('xyz',chargeLocation,'q',q,'R',R);

    for k=1:length(kappa_list)
      kappa = kappa_list(k);


      bemEcfAsym      = makeBemEcfQualMatrices(dielSurfData, pqrData,  epsIn, epsOut);
      bemStern        = makeBemSternMatrices(dielSurfData, sternSurfData, pqrData, ...
				   epsIn, epsOut, kappa);
      phiReacAsym = solveConsistentSternAsym(dielSurfData, sternSurfData, ...
					 pqrData, bemStern, epsIn, epsOut,...
					 kappa, conv_factor, ...
					 asymParams);
      dG_stern(i,j,k) = 0.5 * pqrData.q' * phiReacAsym;
      withstatic_stern(i,j,k) = dG_stern(i,j,k) + sum(pqrData.q)*staticPotential;
      phiReacSym = solveConsistentSternAsym(dielSurfData, sternSurfData, ...
					 pqrData, bemStern, epsIn, epsOut,...
					 kappa, conv_factor, ...
					 symParams);
      dG_Sym_stern(i,j,k) = 0.5 * pqrData.q' * phiReacSym;
      withstatic_sym(i,j,k) = dG_Sym_stern(i,j,k) + sum(pqrData.q)*staticPotential;
    end
  end
end

% from Bardhan12_Jungwirth_Makowski
rscale = 0.92;
sodiumRminOver2 = 1.41075; % new Roux toppar 1.36375;  % standard charmm
sodiumPlus = -93.4 ;

chlorideRminOver2 = 2.27;
chlorideMinus = -95.3 ;

potassiumRminOver2 =1.76375; % new Roux toppar
potassiumPlus  = -73.4;

rubidiumRminOver2 = 1.90;
rubidiumPlus = -66.78 ;

cesiumRminOver2 = 2.1;
cesiumPlus = -60.42 ;

stern_nosalt_anion = squeeze(withstatic_stern(:,1,1));
stern_nosalt_cation = squeeze(withstatic_stern(:,2,1));
stern_nosalt_sym = squeeze(withstatic_sym(:,1,1));

figure; 
plot(R_list, stern_nosalt_anion,'r','linewidth',2);
hold on;
set(gca,'fontsize',16);
xlabel('R_{ion} (Angstrom)');
ylabel('Charging free energy (kcal/mol)');

plot(R_list, stern_nosalt_cation,'b','linewidth',2);
plot(R_list, squeeze(dG_Sym_stern(:,1,1)), 'k--','linewidth',2);

R_list_original = R_list;
load oldnlbc.mat
R_list_old = R_list;
R_list = R_list_original;
plot(R_list_old, L(:,1),'r');
plot(R_list_old, L(:,end),'b');

plot(rscale*[sodiumRminOver2], [sodiumPlus], 'ks','markersize',10,'linewidth',2);

plot(rscale*[potassiumRminOver2], [ potassiumPlus], 'bs','markersize',10,'linewidth',2);

plot(rscale*[rubidiumRminOver2], [rubidiumPlus], 'rs','markersize',10,'linewidth',2);

plot(rscale*[cesiumRminOver2], [cesiumPlus], 'ms','markersize',10,'linewidth',2);

plot(rscale*[chlorideRminOver2], [chlorideMinus], 'k*','markersize',10,'linewidth',2);

legend('New NLBC, q = -1', 'New NLBC, q = +1', 'Poisson', ...
       'Old NLBC, q = -1', 'Old NLBC, q = +1', 'Na, q=+1', ...
       'K, q=+1', 'Rb, q=+1', 'Cs, q=+1', 'Cl, q=-1', ...
		    'Location','EastOutside');

axis([1 2.3 -230 -50])
print -depsc2 withstatic-renormalized.eps
print -dpng withstatic-renormalized.png

