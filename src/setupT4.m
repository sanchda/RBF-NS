%
% Modified set up for the Williamson test case 4.
%
function atm = setupT4(nodes)

% Parameters for the model
atm.sig = 12.74244^2;          % Velocity and geopotential parameter.
atm.g = 9.80616;               % Gravitational constant (m/s^2).
atm.gh0 = 1e5;                 % Mean geopotential (m/s)^2.
atm.a = 6.37122e6;             % Mean radius of the earth (meters).
atm.omega = 7.292e-5;          % Rotation rate of the earth (1/seconds).
atm.u0 = 20;                   % Maximum wind speed.
atm.f0 = 2*atm.omega*sin(pi/4); 
atm.psi0 = -0.03*(atm.gh0/atm.f0);
atm.la0 = 0;                   % longitudinal value for center of initial condition.
atm.th0 = pi/4;                % latitudinal value for center of initial condition.
atm.nu  = 1.789*10^-5;         % dynamic viscosity http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html

atm.pts.x = nodes(:,1); atm.pts.y = nodes(:,2); atm.pts.z = nodes(:,3);
atm.pts.nd = size(atm.pts.x,1);

[atm.pts.la,atm.pts.th,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3)); 
atm.pts.nd = size(atm.pts.x,1);

% Variables for projecting an arbitrary Cartesian vector onto the surface
% of the sphere.
x2 = nodes(:,1).^2; xy = nodes(:,1).*nodes(:,2);
y2 = nodes(:,2).^2; xz = nodes(:,1).*nodes(:,3);
z2 = nodes(:,3).^2; yz = nodes(:,2).*nodes(:,3);
atm.pts.p_u = [1-x2  -xy   -xz];
atm.pts.p_v = [-xy  1-y2   -yz];
atm.pts.p_w = [-xz   -yz  1-z2];

atm.f = 2*atm.omega*sin(atm.pts.th);   % Coriolis force

