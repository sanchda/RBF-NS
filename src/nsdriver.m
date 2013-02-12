cd 'C:\Users\david\Documents\GitHub\RBF-NS\src\';
nu = 1;     % Parameter for the NS equation
eps = 10;   % Shape paramater for the RBF kernel
N = 12;     % Somehow related to the number of centers.  For the ME points,
            % the number of centers is (N+1)^2.
M = 1;      % how many iterations to run the simulation for
h = 0.01;   % timestep

X = getMEPoints(N);
X = X(:,1:3);
% Rotate X through by a small angle
t=0.5;
theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
for i = 1:(N+1)^2
    X(i,:) = (theta*X(i,:)')';
end

% Div-free VFs can be given by curls of scalars.  Just use Qx*grad(u) =
% curl_T(u) on the sphere.

U0 = [-X(:,2) + X(:,3), X(:,1), -X(:,1)];

% u=(x+y+z) => curl = {y - x, z - x, -y + x}
U0 = [X(:,2) - X(:,1), X(:,3) - X(:,1), -X(:,2) + X(:,1)];


onemat = ones((N+1)*(N+1)-1,1);
% Div-only VFs can be given by gradients of scalars.  Px*grad(u):
U0 = [onemat-X(:,1).*X(:,1), -X(:,1).*X(:,2), -X(:,1).*X(:,3)];

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