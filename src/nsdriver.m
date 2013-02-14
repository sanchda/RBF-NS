cd 'C:\Users\david\Documents\GitHub\RBF-NS\src\';
cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';

%==========================================================================
%                         Parameters and Constants                        
%==========================================================================

nu = 1;     % Parameter for the NS equation
omega=1;    % Strength of coriolis force
eps = 10;   % Shape paramater for the RBF kernel
N = 12;     % Somehow related to the number of centers.  For the ME points,
            % the number of centers is (N+1)^2.
N0 = N;     % Highest spherical harmonic in the test
M = 1;      % how many iterations to run the simulation for
h = 0.01;   % timestep

%==========================================================================
%                            Generate RBF nodes                        
%==========================================================================

X = getMEPoints(N);
X = X(:,1:3);
% Rotate X through by a small angle to avoid poles
t=0.5;
theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
for i = 1:(N+1)^2
    X(i,:) = (theta*X(i,:)')';
end

%==========================================================================
%                            Generate initial VF                       
%==========================================================================
% The initial vector field is defined as in GaneshGiaSloan2009 in their
% first numerical trial (eqn 5.1)
%
% U0   = g(0)(W1(x) - W2(x))
% g(t) = nu*exp(-t)(sin(5t) + cos(10t))
% W1   = Z_0^1 + 2Z_1^1
% W2   = Z_0^2 + 2Z_1^2 + 2Z_2^2

Y1 = dsph(1,X(:,1),X(:,2),X(:,3));
Y2 = dsph(2,X(:,1),X(:,2),X(:,3));

U0 = Y1(:,2) + 2*Y1(:,3) - Y2(:,3) - 2*Y2(:,4) - 2*Y2(:,5);
U0 = nu*U0;


%==========================================================================
%                          Generate reference VF                       
%==========================================================================


















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