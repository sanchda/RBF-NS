% AUTHOR:   David Sanchez
% DATE:     August 2012
% MODIFIED: 7/6/2013

% This is the driver script for configuring, initializing, running, and
% visualizing incompressible Navier-Stokes on the sphere, utilizing Grady
% Wright's RBF-based projection method.
%
% TODO:
% * Simplify interface
% * Configure Ganesh's tests not as inlined, but as optional tests
% * Improve interface for visualization (and its initialization)
% * Provide some kind of time-estimate for initialization
% * Improve performance of initialization
%
%cd 'C:\Users\david\Documents\GitHub\RBF-NS\src';

%% Initialization
cd 'C:\Users\david\Desktop\GitHub\RBF-NS\src';

%==========================================================================
%                         Parameters and Constants                        
%==========================================================================

nu = 1/10;      % Parameter for the NS equation
omega = 0;     % Strength of coriolis force
N = 31;        % Somehow related to the number of centers.  For the
               % ME points, the number of centers is (N+1)^2.
N0 = 1;        % Highest spherical harmonic in the test
M = 1;         % how many iterations to run the simulation for
h = 1/(N+1);   % timestep
divFree_geteps = @(N) -0.519226 + 0.106809*(N+1);
eps_Leray = divFree_geteps(N);
beta = 12;
c = (4*pi)^2;
eps_PDE = (beta/c)*(N+1)^(9/8);
surfeps = 2;

disp('Set parameters')
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
X = sortrows(X,3);

disp('RBF nodes loaded')


%==========================================================================
%                            Hessian for matrix RBF                        
%==========================================================================
% I am SO sorry for this.  For ease of use, it was generated using a
% Mathematica script that converts Mathematica output form to Matlab input.
%
% TODO:  generate using Matlab's symbolic toolkit

HGA =  @(x,y,eps) ...
[exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+( ...
-1).*y(3)).^2)).*((-2).*eps.^2+4.*eps.^4.*(x(1)+(-1).*y(1)).^2),4.*exp( ...
1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).* ...
y(3)).^2)).*eps.^4.*(x(1)+(-1).*y(1)).*(x(2)+(-1).*y(2)),4.*exp(1).^((-1).* ...
eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).* ...
eps.^4.*(x(1)+(-1).*y(1)).*(x(3)+(-1).*y(3));4.*exp(1).^((-1).*eps.^2.*(( ...
x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*eps.^4.*(x(1)+( ...
-1).*y(1)).*(x(2)+(-1).*y(2)),exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+( ...
x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*((-2).*eps.^2+4.*eps.^4.*(x(2)+( ...
-1).*y(2)).^2),4.*exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1) ...
.*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*eps.^4.*(x(2)+(-1).*y(2)).*(x(3)+(-1).*y(3)); ...
4.*exp(1).^((-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+( ...
-1).*y(3)).^2)).*eps.^4.*(x(1)+(-1).*y(1)).*(x(3)+(-1).*y(3)),4.*exp(1).^(( ...
-1).*eps.^2.*((x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)) ...
.*eps.^4.*(x(2)+(-1).*y(2)).*(x(3)+(-1).*y(3)),exp(1).^((-1).*eps.^2.*(( ...
x(1)+(-1).*y(1)).^2+(x(2)+(-1).*y(2)).^2+(x(3)+(-1).*y(3)).^2)).*((-2).* ...
eps.^2+4.*eps.^4.*(x(3)+(-1).*y(3)).^2)];

% This is an amazingly wasteful way of defining the operators used
% throughout the NS code, but it has the advantage of providing a
% manifold-agnostic interface.
% TODO: 
% * streamline and clean up
% * determine if any of these have a sparse structure, make sparse
% * vectorize initialization code

[lap projgrad Lx Ly Lz Achol Afull Acrl Adiv PSIfull PSIcrl PSIdiv Pxmat] = nsInitS2(X, HGA, eps_Leray, eps_PDE);
disp('NS working matrices initialized')


%==========================================================================
%                     Optional Test: Vector Laplacian                         
%==========================================================================
%
% I am not proud of this test.
% TODO: swallow regret or fix test
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
U0 = makeGaneshTest1(N0, X, 0, nu);
%U0  = getDivFree(2,X); 
%U0  = U0(:,1:3);

disp('Initial VF generated')


%==========================================================================
%                        Initialize Visualization Stuff                       
%==========================================================================
% 
% Parametrized surface and constants.  -q suffix denotes it's for the
% quiver plot
eps = 1.0;
epsq = 0.1;
phi = @(r2) exp(-eps*eps*r2);
phiq = @(r2) exp(-epsq*epsq*r2);

s2 = @(t,p) [cos(t).*sin(p) sin(t).*sin(p) cos(p)];
sz = [101, 201];
szq = [25 25];
MM = prod(sz);  
MMq = prod(szq);  
NN = size(X,1);

% surface grid parameters
[ll,tt]=meshgrid(linspace(-pi,pi,sz(2)),linspace(-pi,pi,sz(1)));
[llq,ttq]=meshgrid(linspace(-pi,pi,szq(2)),linspace(-pi,pi,szq(1)));
xx=s2(ll(:),tt(:));
xxq=s2(llq(:),ttq(:));

% Interpolate to the grid
re2=(repmat(xx(:,1),[1 NN])-repmat(X(:,1).',[MM 1])).^2;
re2=re2+(repmat(xx(:,2),[1 NN])-repmat(X(:,2).',[MM 1])).^2;
re2=re2+(repmat(xx(:,3),[1 NN])-repmat(X(:,3).',[MM 1])).^2;

re2q=(repmat(xxq(:,1),[1 NN])-repmat(X(:,1).',[MMq 1])).^2;
re2q=re2q+(repmat(xxq(:,2),[1 NN])-repmat(X(:,2).',[MMq 1])).^2;
re2q=re2q+(repmat(xxq(:,3),[1 NN])-repmat(X(:,3).',[MMq 1])).^2;
yy=reshape(xx(:,2),sz);
zz=reshape(xx(:,3),sz);
xx=reshape(xx(:,1),sz);

yyq=reshape(xxq(:,2),szq);
zzq=reshape(xxq(:,3),szq);
xxq=reshape(xxq(:,1),szq);

disp('Visualization parameters set')
%% Simulate with RBF
U = U0;
t=0;

% Initial camera parameters
cpos = [-0.031 -21.392 9.115];
ctarg = [-0.031 0.122 -0.009];
cview = 4.0;

for c = 1:50
% Start with the visualization first to keep from having to handle the
% visualization of the initial condition separately.

% Note that RBFs can't capture constant fields very well, so make sure that
% the field isn't nearly constant (i.e., the zero field) before calling.
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq); 
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq); 
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq); 

% Plot the results
clf
set(gca,'CameraPosition',cpos)
set(gca,'CameraTarget',ctarg)
set(gca,'CameraViewAngle',cview)
%caxis([-0.2,0.2])
colorbar('EastOutside')
cpos = campos;
ctarg = camtarget;
hold on
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
box off
axis off
hold off
text(-1.23,0,1.05,sprintf('t: %f',(c-1)*h),'Fontsize',12)
text(-1.23,0,0.95,sprintf('max: %f',max(max(U))),'Fontsize',12)
text(-1.23,0,0.85,sprintf('min: %f',min(min(U))),'Fontsize',12)
text(-1.23,0,-0.75,sprintf('n: %d',size(U,1)),'Fontsize',12)
text(-1.23,0,-0.85,sprintf('nu: %f',nu),'Fontsize',12)
text(-1.23,0,-0.95,sprintf('omega: %f',omega),'Fontsize',12)

F(c) = getframe(gcf);

[U,t] = navierstokes(X, U, h, t, 1, nu, omega, N0, lap, projgrad, Lx, Ly, Lz, Afull, Acrl, PSIfull, PSIcrl, PSIdiv, Pxmat);

Uganesh = makeGaneshTest1(N0, X, t-h, nu);

err = abs(U - Uganesh);
err = err(:,1).^2 + err(:,2).^2 + err(:,3).^2;
err = sqrt(err);
err = mean(err);

errmat(c) = err;
end


movie(gcf,F,5);
%% Save movie
% TODO: dynamically name these
%
movie2avi(F, 'NS_1024_DAS_7.6.13_trial5_2.avi')

%% View corresponding Ganesh solution
U = U0;

% Initial camera parameters
cpos = [-0.031 -21.392 9.115];
ctarg = [-0.031 0.122 -0.009];
cview = 4.0;

for c = 1:50
t=(c-1)*h;

% Note that RBFs can't capture constant fields very well, so make sure that
% the field isn't nearly constant (i.e., the zero field) before calling.
uu=reshape(phi(re2)*(Achol\(Achol.'\sqrt((U(:,1).^2+U(:,2).^2+U(:,3).^2)))),sz); 
uuq1=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,1))),szq);
uuq2=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,2))),szq);
uuq3=reshape(phiq(re2q)*(Achol\(Achol.'\U(:,3))),szq);

% Plot the results
clf
set(gca,'CameraPosition',cpos)
set(gca,'CameraTarget',ctarg)
set(gca,'CameraViewAngle',cview)
%caxis([-0.2,0.2])
colorbar('EastOutside')
cpos = campos;
ctarg = camtarget;
hold on
quiv = quiver3(xxq,yyq,zzq,uuq1,uuq2,uuq3,1);
htop = surf(xx,yy,zz,uu);
shading interp;
set(htop, 'edgecolor','none')
daspect([1 1 1]);
box off
axis off
hold off
text(-1.23,0,1.05,sprintf('t: %f',t),'Fontsize',12)
text(-1.23,0,0.95,sprintf('max: %f',max(max(U))),'Fontsize',12)
text(-1.23,0,0.85,sprintf('min: %f',min(min(U))),'Fontsize',12)
text(-1.23,0,-0.75,sprintf('n: %d',size(U,1)),'Fontsize',12)
text(-1.23,0,-0.85,sprintf('nu: %f',nu),'Fontsize',12)
text(-1.23,0,-0.95,sprintf('omega: %f',omega),'Fontsize',12)

G(c) = getframe(gcf);

t = t+h;
U = makeGaneshTest1(N0, X, t, nu);
end

movie(gcf,G,1)

%% Save movie
% TODO: dynamically name these
%
movie2avi(G, 'NS_1024_GAN_7.6.13_trial5.avi')
