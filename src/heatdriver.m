cd 'C:\Users\david\Documents\GitHub\RBF-NS\src\';
alpha = 1;  % Parameter for the heat equation
eps = 1;    % Shape paramater for the RBF kernel
N = 5;     % Somehow related to the number of centers.  For the ME points,
            % the number of centers is (N+1)^2.
M = 1;      % how many iterations to run the simulation for
h = 0.01;   % timestep

X = getMEPoints(N);
X = X(2:end,:);  % This first point causes problems for some reason.

%DEBUG:  if one makes some shitty argument about the point being approached
%to z->1 from (x,y)=(0,0), the limit works out for ej and dj.  The reason
%why the original point breaks is because this (and the antipodal point,
%I guess) is the only point where the denominator goes to zero.

% Div-free VFs can be given by curls of scalars.  Just use Qx*grad(u) =
% curl_T(u) on the sphere.

% u=x => curl = {0, z, -y}
U0 = [0.*X(:,1), X(:,3), -X(:,2)];

% u=(x+y+z) => curl = {y - x, z - x, -y + x}
U0 = [X(:,2) - X(:,1), X(:,3) - X(:,1), -X(:,2) + X(:,1)];


% Debug!
% heateqn should be rewritten to remove the assumption that the points in X
% lie on the sphere.  This test is to see how the differential operators
% act on the plane.
% N=600;
% p = haltonset(2);
% p = haltonset(3,'Skip',1e3,'Leap',1e2);
% x = net(p,N);
% z = zeros(N,1);
% X = [x(:,1) x(:,2) z];
% U0 = X(:,1).*X(:,1);
% 

U = heateqn(U0, X, h, M, alpha, eps);
% 
% title('function of x');
% subplot(2,1,1);
% scatter3(X(:,1),X(:,2),U0,30,U0,'.')
% 
% subplot(2,1,2);
% scatter3(X(:,1),X(:,2),U,30,U,'.')
% title('derivative');


% Convert the Cartesian coordinates to longitudinal/latitudinal coordinates
% for visualization
[p,t] = cart2sph(X(:,1), X(:,2), X(:,3));

% Display U according to p,t coordinates
subplot(2,3,1);
scatter3(p,t,U(:,1),30,U(:,1),'.')
subplot(2,3,2);
scatter3(p,t,U(:,2),30,U(:,2),'.')
subplot(2,3,3);
scatter3(p,t,U(:,3),30,U(:,3),'.')

subplot(2,3,4);
scatter3(p,t,U0(:,1),30,U0(:,1),'.')
subplot(2,3,5);
scatter3(p,t,U0(:,2),30,U0(:,2),'.')
subplot(2,3,6);
scatter3(p,t,U0(:,3),30,U0(:,3),'.')