cd 'C:\Users\david\Documents\GitHub\RBF-heat\src\';
alpha = 1;  % Parameter for the heat equation
eps = 0.5;  % Shape paramater for the RBF kernel
N = 10;     % Somehow related to the number of centers.  For the ME points,
            % the number of centers is (N+1)^2.
M = 1;      % how many iterations to run the simulation for
h = 0.001;  % timestep

X = getMEPoints(N);

% Generate some initial data.
U0 = zeros(size(X,1),1);
U0 = U0 + 1;

U = heateqn(U0, X, h, M, alpha, eps);

% Convert the Cartesian coordinates to longitudinal/latitudinal coordinates
% for visualization
[p,t] = cart2sph(X(:,1), X(:,2), X(:,3));

% Display U according to p,t coordinates
scatter3(p,t,U,30,U,'.')
