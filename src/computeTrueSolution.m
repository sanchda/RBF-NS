% True solution to the modified test case 4.  uc contains the velocity in
% Cartesian coordinates with uc(:,1:3) the x,y,z direction, respectively.
function uc = computeTrueSolution(atm,t)

% Extract out the constants from the atm structure to make the code easier
% to read.
sig = atm.sig;
a = atm.a;
u0 = atm.u0;
psi0 = atm.psi0;
la0 = atm.la0;
th0 = atm.th0;
x = atm.pts.x;
y = atm.pts.y;
z = atm.pts.z;

% Compute the center of the low pressure system at time t.
x0 = cos(la0+u0/a*t).*cos(th0); y0 = sin(la0+u0/a*t).*cos(th0); z0 = sin(th0);

% Dot product between nodes and the center of the low pressure.  This is 
% equation (129) from Williamson expressed in Cartesian coordinates.
dp = x*x0 + y*y0 + z*z0; 

% Cross product between nodes and the center of the low pressure.
cp = [y*z0-z*y0 z*x0-x*z0 x*y0-y*x0];
b = sig*(1-dp)./(1+dp);
psibar = psi0*exp(-b);
tempu = 2*sig*psibar./(a*(1 + dp).^2);

% Flow field generated from the low pressure (eq (123) from Williamson without
% the ubar term).
uc = cp.*repmat(tempu,[1 3]);
tempu = u0*(2*z).^14.*(1-z.^2).^(13/2);

% Baseline flow field (equation (126) from Williamson, but in Cartesian coords.)
ub(:,1) = -tempu.*y;
ub(:,2) = tempu.*x;
ub(:,3) = 0;

% Full flow field (flow from low pressure, along with the baseline).
uc(:,1) = uc(:,1) + ub(:,1);
uc(:,2) = uc(:,2) + ub(:,2);

end
