function [A, P] = calcYoonSternAsymPert(asymParams, dielSurfData, ...
					sternSurfData, bem, pqr, ...
					epsIn, epsOut, kappa, phiDielBndy, ...
					dphiDnDielBndy, phiSternBndy, ...
					dphiDnSternBndy)

alpha = asymParams.alpha;
beta  = asymParams.beta;
EfieldOffset = asymParams.EfieldOffset;
deltaOffset  = -alpha * tanh(beta*0-EfieldOffset);

Efield = -bem.dielChargeOp.dphidnCoul * pqr.q - bem.dielDielOp.K'*dphiDnDielBndy ...
	 + bem.dielDielOp.W*phiDielBndy;

h = (alpha*(tanh(beta*Efield-EfieldOffset)) +deltaOffset);
f = (epsIn/(epsOut-epsIn)) - h;

r1 = -sum(pqr.q)/epsOut;
r2 = sum(sternSurfData.weights .* dphiDnSternBndy);
r1Overr2 = r1/r2; %epsOut/epsOut;

A = [bem.A11 bem.A12 bem.A13 bem.A14;
     bem.A21 bem.A22_base*diag(f./(1+f)) bem.A23 bem.A24;
     bem.A31 bem.A32_base*diag(f./(1+f)) bem.A33 bem.A34;
     bem.A41 bem.A42 bem.A43 bem.A44_base*(r1Overr2)];

P = sparse([diag(diag(bem.A11)) diag(diag(bem.A12)) 0*bem.A13 0*bem.A14;
     diag(diag(bem.A21)) diag(diag(bem.A22_base))*diag(f./(1+f)) ...
     0*bem.A23 0*bem.A24;
     0*bem.A31 0*bem.A32_base diag(diag(bem.A33)) diag(diag(bem.A34));
     0*bem.A41 0*bem.A42 diag(diag(bem.A43)) diag(diag(bem.A44_base))*(r1Overr2)]);
