% Simplify the execution by changing current directory to the script's
% local directory.  There's probably a better way, but for now this should
% be sufficient (given one uses a sufficiently nice pathname).
[cur_dir,~,~] = fileparts(which('nsdriver.m'));
cd(cur_dir);

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

%==============================debug parameters============================
% Check the version of the vector Laplacian saved in testVecLap.m
test_lap = 0;

%==========================================================================
%                            Generate RBF nodes                        
%==========================================================================

X = getMEPoints(N);
W = X(:,4);
X = X(:,1:3);
% Rotate X through by a small angle to avoid poles
t=0.5;
theta = [1 0 0;0 cos(t) -sin(t);0 sin(t) cos(t)];
for i = 1:(N+1)^2
    X(i,:) = (theta*X(i,:)')';
end

%==========================================================================
%                     Optional Test: Vector Laplacian                         
%==========================================================================
if test_lap == 1
    
   divFree_L = 1;                                     % how high is the degree of the VSH
   divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);  % empirically fitted shape parameter for the vector Laplacian
   divFree_U0  = getDivFree(2,X); 
   divFree_U0  = divFree_U0(:,1:3);                   % use a divergence-free VSH
   divFree_eps = divFree_geteps(N);
   [~, divFree_Maxerr, divFree_L2err, divFree_kappa] = testVecLap(X, W, divFree_U0, divFree_eps);
   [divFree_L2err divFree_Maxerr]
   
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

g  = @(t) nu*exp(-t)*(sin(5*t)+cos(10*t));

Z = getDivFree(1,X);
W1 = [Z(:,4) Z(:,5) Z(:,6)] + 2*[Z(:,7) Z(:,8) Z(:,9)];

Z = getDivFree(2,X);
W2 = [Z(:,7) Z(:,8) Z(:,9)] + 2*[Z(:,10) Z(:,11) Z(:,12)] + 2*[Z(:,13) Z(:,14) Z(:,15)];

clear(W1,W2,Z);
U = nu*(W1 - W2);


%==========================================================================
%                          Timestep + Check                       
%==========================================================================
% 
% Run the simulation and check the output against the reference timestep in
% makeGaneshTest1.  Note that this test has no explicit Coriolis force.
%

U = navierstokes(X,U,H,h,1,epsilon,nu,0);
U_ref = makeGaneshTest1(N0, X, c*h, nu);










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