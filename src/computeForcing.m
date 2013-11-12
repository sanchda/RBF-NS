% Computes the forcing terms for the modified Williamson's Test Case 4. All
% computations are for the Cartesian form of the Navier-Stokes equations.
% F(:,1) = x forcing, F(:,2) = y forcing, F(:,3) = z forcing.
function F = computeForcing(atm,t)

uc = computeTrueSolution(atm,t);

% Extract out the constants from the atm structure to make the code easier
% to read.
sig = atm.sig;
a = atm.a;
u0 = atm.u0;
psi0 = atm.psi0;
la0 = atm.la0;
th0 = atm.th0;
f = atm.f;
x = atm.pts.x;
y = atm.pts.y;
z = atm.pts.z;

z2 = z.^2;     % Square z and store it since we use it in several places.

% Compute the center of the low pressure system at time t.
x0 = cos(la0+u0/a*t).*cos(th0); y0 = sin(la0+u0/a*t).*cos(th0); z0 = sin(th0);

% Compute the cross product of x and xc
cp = [y*z0-z*y0 z*x0-x*z0 x*y0-y*x0];
cpdt = u0*[-x0*z -y0*z (x*x0 + y*y0)]/a;

% Compute the dot product of x and xc
dp = x*x0 + y*y0 + z*z0;
dpp = 1 + dp;
dpp2 = (1 + dp).^2;
dpp4 = (1 + dp).^4;
dpm = 1-dp;

% Compute psibar
b = sig*(dpm./dpp);
psibar = psi0*exp(-b);
dpsibar = 2*sig*psibar./(a*dpp2);
d2psibar = 4*(sig-dpp)*sig.*psibar./(a^2*dpp4);

%
% Compute the time derivatives.
%

% Momentum equations.
F = cpdt.*repmat(dpsibar,[1 3]) + cp.*repmat(u0*(y*x0-x*y0).*d2psibar,[1 3]);

%
% Compute the u.grad(u) terms in the momentum equations.
%

% Temporary variables:
temp = d2psibar.*(x0*uc(:,1) + y0*uc(:,2) + z0*uc(:,3));
temp2 = u0/a*(2*z).^14.*(1-z2).^(13/2);
dtemp2 = u0/a*(28*(2*z).^13.*(1-z2).^(13/2) - 13*(2*z).^14.*(1-z2).^(11/2).*z);

F(:,1) = F(:,1) + cp(:,1).*temp + z0/a*dpsibar.*uc(:,2) - y0/a*dpsibar.*uc(:,3) - temp2.*uc(:,2) - dtemp2.*y.*uc(:,3);
F(:,2) = F(:,2) + cp(:,2).*temp - z0/a*dpsibar.*uc(:,1) + x0/a*dpsibar.*uc(:,3) + temp2.*uc(:,1) + dtemp2.*x.*uc(:,3);
F(:,3) = F(:,3) + cp(:,3).*temp + y0/a*dpsibar.*uc(:,1) - x0/a*dpsibar.*uc(:,2);

% 
% Compute the Lagrange multiplier term (for constraining the flow to the surface)
% of the sphere.
%
nrmu = uc(:,1).^2 + uc(:,2).^2 + uc(:,3).^2;
F(:,1) = F(:,1) + 1/a*x.*nrmu;
F(:,2) = F(:,2) + 1/a*y.*nrmu;
F(:,3) = F(:,3) + 1/a*z.*nrmu;

% 
% Compute the Coriolis force in the momemtum equations.
%
F(:,1) = F(:,1) + f.*(-z.*uc(:,2) + y.*uc(:,3));
F(:,2) = F(:,2) + f.*( z.*uc(:,1) - x.*uc(:,3));
F(:,3) = F(:,3) + f.*(-y.*uc(:,1) + x.*uc(:,2));

end

